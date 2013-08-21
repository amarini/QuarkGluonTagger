

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

#define TREEMC "/afs/cern.ch/work/p/pandolf/public/sunilFlat_DiJet_flatQCD_P6_Dijets_12Aug_ptHatWeight.root"
#define TREEDATA "/afs/cern.ch/work/p/pandolf/public/sunilFlat_DiJet_data2012ABCD_MBPD_12Jul.root"

#include "BaseDoubleMin.C"
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

//const inline bool operator<(TRIO&x,TRIO&y){return x.value<y.value;}

class Analyzer:public BaseAnalyzer{
public:
	Analyzer(){
		qgl=new QGLikelihoodCalculator("/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/data/");//ReducedHisto_2012.root");
		qgmlp=new QGMLPCalculator("MLP","/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/data/",true); //prob
		}
	void Loop(TChain *t,int type); //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc
	void LoadBins();
//private:
	//histo reweight
//--- QGL QGMLP
	QGLikelihoodCalculator *qgl;
	QGMLPCalculator *qgmlp;
};


void Analyzer::Loop(TChain *t,int type){ //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc

		
	//	treeVarInt["nvtx"] = -999; t->SetBranchAddress("nvtx",&treeVarInt["nvtx"]);
		treeVar["rhoPF"] = -999; t->SetBranchAddress("rhoPF",&treeVar["rhoPF"]);
		treeVar["rho"] = -999; t->SetBranchAddress("rhoPF_allJets",&treeVar["rho"]);
		Float_t jetPt[4]; t->SetBranchAddress("ptJet",&jetPt);
		Float_t jetEta[4]; t->SetBranchAddress("etaJet",&jetEta);
		//Float_t jetAxis1_QC[4]; t->SetBranchAddress("axis1_QCJet",&jetAxis1_QC);
		//Float_t jetAxis2_QC[4]; t->SetBranchAddress("axis2_QCJet",&jetAxis2_QC);
		//Float_t jetPtD_QC[4]; t->SetBranchAddress("ptD_QCJet",&jetPtD_QC);
		//Int_t jetChgPart_QC[4]; t->SetBranchAddress("nChg_QCJet",&jetChgPart_QC);
		//Int_t jetNeutralPart_ptcut[4]; t->SetBranchAddress("nNeutral_ptcutJet",&jetNeutralPart_ptcut);
		Int_t pdgIdJet[4]; t->SetBranchAddress("pdgIdJet",&pdgIdJet);
		treeVar["eventWeight"] = -999; t->SetBranchAddress("eventWeight",&treeVar["eventWeight"]);
		Float_t jetQGL[4]; t->SetBranchAddress("qglJet",&jetQGL);
		Float_t jetMLP[4]; t->SetBranchAddress("qgMLPJet",&jetMLP);
		
		
		if(type&4) {lmin=1.0;lmax=0;} //reset lmin-lmax
		if(type&1) {delete h_data; CreateHisto(1);varAllData.clear();}
		if(type&10) {delete h_mc; CreateHisto(2);} //8+2
		if(type&32) {varAll.clear();} //reset varAll

		for(int i=0;i<t->GetEntries() ;i++)
			{
			t->GetEntry(i);
			treeVar["ptJet0"]=jetPt[0];
			treeVar["ptJet1"]=jetPt[1];
			treeVar["etaJet0"]=jetEta[0];
			treeVar["etaJet1"]=jetEta[1];
			treeVar["rhoPFJets"]=treeVar["rho"];
		
			bool sel1=true;
			bool sel0=true;
			//TP Pt 1<->0
			if((treeVar["ptJet1"]<PtMin)||(treeVar["ptJet1"]>PtMax)||(fabs(treeVar["etaJet0"])<EtaMin)||(fabs(treeVar["etaJet0"])>EtaMax)|| (treeVar["rhoPF"]<RhoMin)||(treeVar["rhoPF"]>RhoMax))sel0=false;
			if((treeVar["ptJet0"]<PtMin)||(treeVar["ptJet0"]>PtMax)||(fabs(treeVar["etaJet1"])<EtaMin)||(fabs(treeVar["etaJet1"])>EtaMax)|| (treeVar["rhoPF"]<RhoMin)||(treeVar["rhoPF"]>RhoMax))sel1=false;
			
			if(sel1==false && sel0==false) continue;
			
			treeVarInt["pdgIdPartJet0"]=pdgIdJet[0];
			treeVarInt["pdgIdPartJet1"]=pdgIdJet[1];
				
			treeVar["QGLMLP"]=jetMLP[0];
			treeVar["QGLHisto"]=jetQGL[0];
			{
		//	treeVar["eventWeight"]=1.;
			treeVar["PUReWeight"]=1.; //already contained in eventW
			}
		
			if(type&1){
				string var=varName;
				alpha=1; beta=0;
					treeVar["eventWeight"]=1.; //data
					treeVar["PUReWeight"]=1;
				if(sel0){
				treeVar["QGLMLP"]=jetMLP[0];
				treeVar["QGLHisto"]=jetQGL[0];
				FillHisto(h_data,var);
				varAllData.push_back(treeVar[var]);
				}
				if(sel1){
				treeVar["QGLMLP"]=jetMLP[1];
				treeVar["QGLHisto"]=jetQGL[1];
				FillHisto(h_data,var);
				varAllData.push_back(treeVar[var]);
				}
				}	
			if(type&2){
				//mc
				if(sel0){
				treeVar["QGLMLP"]=jetMLP[0];
				treeVar["QGLHisto"]=jetQGL[0];
				FillHisto(h_mc,varName);
				}
				if(sel1){
				treeVar["QGLMLP"]=jetMLP[1];
				treeVar["QGLHisto"]=jetQGL[1];
				FillHisto(h_mc,varName);
				}
				}
			if(type&4){
				string var=varName;
				if(sel1){
				treeVar["QGLMLP"]=jetMLP[1];
				treeVar["QGLHisto"]=jetQGL[1];
				if(lmin>treeVar[var]) lmin=treeVar[var];
				if(lmax<treeVar[var]) lmax=treeVar[var];
				}if(sel0){
				treeVar["QGLMLP"]=jetMLP[0];
				treeVar["QGLHisto"]=jetQGL[0];
				if(lmin>treeVar[var]) lmin=treeVar[var];
				if(lmax<treeVar[var]) lmax=treeVar[var];
				}
				}
			if(type&8){
				alpha=1;beta=0;
				if(sel0){
				treeVar["QGLMLP"]=jetMLP[0];
				treeVar["QGLHisto"]=jetQGL[0];
				if( treeVarInt["pdgIdPartJet0"] ==21) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet0"]) < 5) {alpha=a_q;beta=b_q;}
				if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet0"])== 999) {alpha=a_g;beta=b_g;}
				FillHisto(h_mc,varName);
				}
				if(sel1){
				treeVar["QGLMLP"]=jetMLP[1];
				treeVar["QGLHisto"]=jetQGL[1];
				if( treeVarInt["pdgIdPartJet1"] ==21) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet1"]) < 5) {alpha=a_q;beta=b_q;}
				if( fabs(treeVarInt["pdgIdPartJet1"])== 0) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet1"])== 999) {alpha=a_g;beta=b_g;}
				FillHisto(h_mc,varName);
				}
				}
			if(type&16){
				for(int j=0;j<int(alphaFast.size());j++)
					{
					alpha=alphaFast[j];
					beta=betaFast[j];
					if(sel0){
					treeVar["QGLMLP"]=jetMLP[0];
					treeVar["QGLHisto"]=jetQGL[0];
					FillHisto(h_mcFast[j],varName);
					}
					if(sel1){
					treeVar["QGLMLP"]=jetMLP[1];
					treeVar["QGLHisto"]=jetQGL[1];
					FillHisto(h_mcFast[j],varName);
					}
					}
				}
			if(type&32){ 
					float weight=treeVar["eventWeight"]*treeVar["PUReWeight"];
					if(sel0){
					treeVar["QGLMLP"]=jetMLP[0];
					treeVar["QGLHisto"]=jetQGL[0];
					varAll.push_back(TRIO(treeVarInt["pdgIdPartJet0"],treeVar[varName],weight)); //w=-1 default
					}
					if(sel1){
					treeVar["QGLMLP"]=jetMLP[1];
					treeVar["QGLHisto"]=jetQGL[1];
					varAll.push_back(TRIO(treeVarInt["pdgIdPartJet1"],treeVar[varName],weight)); //w=-1 default
					}
				}
			}
}


void Analyzer::LoadBins(){
		ResetBins();
		PtBins.push_back(  pair<float,float>(30,40) );
		PtBins.push_back(  pair<float,float>(40,50) );
		PtBins.push_back(  pair<float,float>(50,65) );
		PtBins.push_back(  pair<float,float>(65,85) );
		PtBins.push_back(  pair<float,float>(85,110) );
		PtBins.push_back(  pair<float,float>(110,140) );
		PtBins.push_back(  pair<float,float>(140,180) );
		PtBins.push_back(  pair<float,float>(180,230) );
		PtBins.push_back(  pair<float,float>(230,300) );
		PtBins.push_back(  pair<float,float>(300,4000) );
		
		RhoBins.push_back(  pair<float,float>(0,100) );
		
		EtaBins.push_back(  pair<float,float>(0,2) );
	//	EtaBins.push_back(  pair<float,float>(2.0,2.5) );
	//	EtaBins.push_back(  pair<float,float>(2.5,3.0) );
		EtaBins.push_back(  pair<float,float>(3,4.7) );
	}



//-----------------------------------DOUBLE MIN FAST ---------------------------------------------------------


//----------------------------------------------------------------------------------------------------



int ComputeDoubleMinDiJetTP(){
	system("[ -f output.root ] && rm output.root");
	TChain *mc=new TChain("tree_passedEvents");
	TChain *data=new TChain("tree_passedEvents");
		mc->Add(TREEMC);

		data  ->Add(TREEDATA);
	Analyzer A;
	A.nstep=20;
	A.compress=1;
	A.varName="QGLHisto";
//	A.varName="QGLMLP";
	A.CreateHisto();
	A.SetTrees(mc,data);
	//fprintf(stderr,"Start Points from ZJet\n");
	A.readStartPoints=0;
	//A.filenameStartPoints="../data/SystZJetHbb_2013_07_23.txt" ;
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
	//A.SpanMin();
	return 0;
	}

Analyzer *Check(float PtMin,float PtMax,float RhoMin,float RhoMax, float EtaMin,float EtaMax, float alpha, float beta , const char * varName="QGLHisto",float lmin=0,float lmax=1,TCanvas **pC=NULL){
Analyzer *A=new Analyzer();
TChain *mc=new TChain("tree_passedEvents");
TChain *data=new TChain("tree_passedEvents");
mc->Add(TREEMC);
data  ->Add(TREEDATA);
A->nstep=15; A->varName=varName;
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
	hmc1->GetYaxis()->SetRangeUser(0,0.30);

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

Analyzer *CheckDouble(float PtMin,float PtMax,float RhoMin,float RhoMax, float EtaMin,float EtaMax, float a_q, float b_q, float a_g, float b_g , const char * varName="QGLHisto",float lmin=0,float lmax=1,TCanvas **pC=NULL){
Analyzer *A=new Analyzer();
TChain *mc=new TChain("tree_passedEvents");
TChain *data=new TChain("tree_passedEvents");
cout<<"Opening Files "<<TREEMC<<" "<<TREEDATA<<endl;
mc->Add(TREEMC);
data  ->Add(TREEDATA);
A->nstep=15; A->varName=varName;
A->RhoMin=RhoMin; A->RhoMax=RhoMax;A->PtMin=PtMin;A->PtMax=PtMax; A->EtaMin=EtaMin;A->EtaMax=EtaMax;
A->CreateHisto();
A->SetTrees(mc,data);
A->alpha=1;
A->beta=0;
A->Loop(mc,2);
	A->compress=1;
A->Loop(data,1);
TH1F* h_mc0=(TH1F*)A->h_mc->Clone("h_mc0");h_mc0->SetLineColor(kGreen);
A->lmin=lmin;
A->lmax=lmax;
//A->alpha=alpha;
//A->beta=beta;
//A->Loop(mc,2);
//DoubleMin
A->Loop(mc,32);
A->a_q=a_q;A->b_q=b_q;A->a_g=a_g;A->b_g=b_g;
A->LoopFast();

gStyle->SetOptStat(0);
TCanvas *c=new TCanvas("c","c",800,800);
A->h_mc->SetLineColor(kBlue);
	A->h_mc->GetXaxis()->SetTitle(varName);
	A->h_mc->SetMinimum(0);
TH1F * hmc1=(TH1F*)A->h_mc->DrawNormalized("HIST");
	hmc1->GetYaxis()->SetRangeUser(0,0.30);

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

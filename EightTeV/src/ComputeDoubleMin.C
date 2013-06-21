

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
#include "BaseDoubleMin.C"

using namespace std;

class Analyzer: public BaseAnalyzer {
public:
	Analyzer(): BaseAnalyzer() {}
	void Loop(TChain *t,int type); //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc
	void LoadBins();
};

void Analyzer::Loop(TChain *t,int type){ //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc

		treeVar["betaStarJet0"]		=-999;t->SetBranchAddress("betaStarJet0",	&treeVar["betaStarJet0"]	);
		treeVarInt["nvertex"]		=-999;t->SetBranchAddress("nvertex",		&treeVarInt["nvertex"]		);
		treeVar["deltaPhi_jet"]		=-999;t->SetBranchAddress("deltaPhi_jet",	&treeVar["deltaPhi_jet"]	);
		treeVar["mZ"]			=-999;t->SetBranchAddress("mZ",			&treeVar["mZ"]			);
		treeVarInt["nPFCand_QC_ptCutJet"]=-999;t->SetBranchAddress("nPFCand_QC_ptCutJet",&treeVarInt["nPFCand_QC_ptCutJet"]			);
		treeVar["ptD_QCJet0"]			=-999;t->SetBranchAddress("ptD_QCJet0",			&treeVar["ptD_QCJet0"]			);
		treeVar["axis1_QCJet0"]			=-999;t->SetBranchAddress("axis1_QCJet0",			&treeVar["axis1_QCJet0"]			);
		treeVar["axis2_QCJet0"]			=-999;t->SetBranchAddress("axis2_QCJet0",			&treeVar["axis2_QCJet0"]			);
		treeVar["ptZ"]			=-999;t->SetBranchAddress("ptZ",		&treeVar["ptZ"]			);
		treeVar["ptJet0"]		=-999;t->SetBranchAddress("ptJet0",		&treeVar["ptJet0"]		);
		treeVar["rhoPF"]		=-999;t->SetBranchAddress("rhoPF",		&treeVar["rhoPF"]		);
		treeVar["etaJet0"]		=-999;t->SetBranchAddress("etaJet0",		&treeVar["etaJet0"]		);
		treeVarInt["pdgIdPartJet0"]	=-999;if(type!=1)t->SetBranchAddress("pdgIdPartJet0",	&treeVarInt["pdgIdPartJet0"]); //only mc 

		treeVar["QGLHisto"]		=-999;t->SetBranchAddress("QGLHisto",		&treeVar["QGLHisto"]		);
		treeVar["QGLMLP"]		=-999;t->SetBranchAddress("QGLMLP",		&treeVar["QGLMLP"]		);
		treeVar["QGLHistoFwd"]		=-999;if(type==1)t->SetBranchAddress("QGLHistoFwd",	&treeVar["QGLHistoFwd"]		); //only data
		treeVar["QGLMLPFwd"]		=-999;if(type==1)t->SetBranchAddress("QGLMLPFwd",	&treeVar["QGLMLPFwd"]		); //only data
		
		treeVar["PUReWeight"]		=1   ;if(type!=1)t->SetBranchAddress("PUReWeight",	&treeVar["PUReWeight"]		); //only mc
		treeVar["eventWeight"]		=1   ;if(type!=1)t->SetBranchAddress("eventWeight",	&treeVar["eventWeight"]		); //only mc
		
		if(type&4) {lmin=1.0;lmax=0;} //reset lmin-lmax
		if(type&1) {delete h_data; CreateHisto(1);}
		if(type&10) {delete h_mc; CreateHisto(2);} //8+2
		if(type&32) {varAll.clear();} //reset varAll

		for(int i=0;i<t->GetEntries();i++) 
			{
			t->GetEntry(i);
			//printf("Count -1: %f %f %f\n",treeVar["ptJet0"],treeVar["etaJet0"],treeVar["rhoPF"]);
			if((treeVar["ptJet0"]<PtMin)||(treeVar["ptJet0"]>PtMax)||(fabs(treeVar["etaJet0"])<EtaMin)||(fabs(treeVar["etaJet0"])>EtaMax)|| (treeVar["rhoPF"]<RhoMin)||(treeVar["rhoPF"]>RhoMax))continue;
			//printf("Count 0\n");
			//Z
			if( (treeVar["mZ"]<70) || (treeVar["mZ"]>110)|| ( fabs(treeVar["deltaPhi_jet"])<3.1415-0.5) ) continue;
			//printf("Count 1 -- Z\n");
			if( (EtaMin<2.5) && (treeVar["betaStarJet0"] > 0.2 * TMath::Log( treeVarInt["nvertex"] - 0.67))  ) continue;
			//printf("Count 2 -- beta*\n");
			if(( (treeVar["axis1_QCJet0"] <=0 ) || (treeVar["axis2_QCJet0"] <=0) ))continue;
			//printf("Count 3 -- axis\n");
			if( treeVarInt["nPFCand_QC_ptCutJet"] <=0 )continue;
			if( treeVar["ptD_QCJet0"] <=0 )continue;
			//printf("Count 4 -- mult\n");
			//---------------------------
		
			if(type&1){
				//printf("passed selection - type 1 --> %.3f - %.3f\n",treeVar[varName],treeVar[varName+"Fwd"]);
				string var=varName;
				if(EtaMin>2.5)var+="Fwd"; //only in data fwd
				alpha=1; beta=0;
				FillHisto(h_data,var);
				}	
			if(type&2){
				//mc
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
					//varAll.push_back(pair<int,float>(treeVarInt["pdgIdPartJet0"],treeVar[varName]));
					float weight=treeVar["eventWeight"]*treeVar["PUReWeight"];
					varAll.push_back(TRIO(treeVarInt["pdgIdPartJet0"],treeVar[varName],weight)); //w=-1 default
				}
			}
}

void Analyzer::LoadBins(){
		PtBins.push_back(  pair<float,float>(30,50) );
		PtBins.push_back(  pair<float,float>(50,80) );
		PtBins.push_back(  pair<float,float>(80,120) );
		PtBins.push_back(  pair<float,float>(120,250) );
		
		RhoBins.push_back(  pair<float,float>(0,15) );
		RhoBins.push_back(  pair<float,float>(15,40) );
		
		EtaBins.push_back(  pair<float,float>(0,2) );
		EtaBins.push_back(  pair<float,float>(3,4.7) );
	}


int ComputeDoubleMin(){
	system("[ -f output.root ] && rm output.root");
	TChain *mc=new TChain("tree_passedEvents");
	TChain *data=new TChain("tree_passedEvents");
		data->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DoubleMu*.root");
		data->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DoubleE*.root");

		mc  ->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DYJetsToLL_M-50*.root");
	Analyzer A;
	A.stp0=.002;A.stp1=0.002;
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

#include "TStyle.h"
#include "TLegend.h"
Analyzer *Check(float PtMin,float PtMax,float RhoMin,float RhoMax, float EtaMin,float EtaMax, float alpha, float beta , const char * varName="QGLHisto",float lmin=0,float lmax=1,TCanvas **pC=NULL){
Analyzer *A=new Analyzer();
	TChain *mc=new TChain("tree_passedEvents");
	TChain *data=new TChain("tree_passedEvents");
		data->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DoubleMu*.root");
		data->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DoubleE*.root");

		mc  ->Add("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DYJetsToLL_M-50*.root");
A->nstep=20; A->varName=varName;
A->RhoMin=RhoMin; A->RhoMax=RhoMax;A->PtMin=PtMin;A->PtMax=PtMax; A->EtaMin=EtaMin;A->EtaMax=EtaMax;
A->CreateHisto();
A->SetTrees(mc,data);
A->alpha=1;
A->beta=0;
A->Loop(mc,2);
A->Loop(data,1);
TH1F* h_mc0=(TH1F*)A->h_mc->Clone("h_mc0");h_mc0->SetLineColor(kGreen);
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



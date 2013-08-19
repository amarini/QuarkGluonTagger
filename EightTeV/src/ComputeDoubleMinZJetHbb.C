

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
#include "TLorentzVector.h"

//---compute QGL & QGLMPL on the fly
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/QGMLPCalculator.cc"
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/parameters.cc"
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/QGLikelihoodCalculator.cc"
#include "/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/src/Bins.cc"

#include "BaseDoubleMin.C"

//#define DATATREE "/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012ABCD_17Jan.root"
//#define MCTREE "/afs/cern.ch/work/s/sunil/public/forTom/analysis_Zjets_17Jan.root"

#define DATATREE "/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012ABCD_MuPD_12Jul.root"
#define MCTREE "/afs/cern.ch/work/s/sunil/public/forTom/analysis_Zjets_12Jul.root"
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
		TFile *Fpuw=TFile::Open("/afs/cern.ch/work/s/sunil/public/forTom/PUhist.root");
//		TFile *Fptetaw=TFile::Open("/afs/cern.ch/work/s/sunil/public/forTom/Jetpteta_rewt2D_flatP6.root ");
		puw=(TH1F*)Fpuw->Get("hist_puWT")->Clone("hPU_wt");
	//	ptetaw=(TH2F*)Fptetaw->Get("hist_WT")->Clone("hPtEta_wt");
		qgl=new QGLikelihoodCalculator("/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/data/");//ReducedHisto_2012.root");
		qgmlp=new QGMLPCalculator("MLP","/afs/cern.ch/user/a/amarini/work/CMSSW_5_3_6/src/QuarkGluonTagger/EightTeV/data/",true); //prob
		}
	void Loop(TChain *t,int type); //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc
	void LoadBins();
//private:
	//histo reweight
		TH1F* puw;
		TH2F* ptetaw;
//--- QGL QGMLP
	QGLikelihoodCalculator *qgl;
	QGMLPCalculator *qgmlp;
};


void Analyzer::Loop(TChain *t,int type){ //type|=4 : compute lmin,lmax; type|=1 data type |=2 mc

		
		treeVarInt["nvtx"] = -999; t->SetBranchAddress("nvtx",&treeVarInt["nvtx"]);
		treeVar["rho"] = -999; t->SetBranchAddress("rho",&treeVar["rho"]);
		treeVar["rhoIso"] = -999; t->SetBranchAddress("rhoIso",&treeVar["rhoIso"]);
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

		Float_t         MuPt[2];   
  		Float_t         MuEta[2];
  	 	Float_t         MuPhi[2];
   		Float_t         MuEnergy[2];
   		Float_t         Invmass;
		   t->SetBranchAddress("MuPt", &MuPt);
		   t->SetBranchAddress("MuEta", &MuEta);
		   t->SetBranchAddress("MuPhi", &MuPhi);
		   t->SetBranchAddress("MuEnergy", &MuEnergy);
		   t->SetBranchAddress("Invmass", &Invmass);
		
		
		if(type&4) {lmin=1.0;lmax=0;} //reset lmin-lmax
		if(type&1) {delete h_data; CreateHisto(1);varAllData.clear();}
		if(type&10) {delete h_mc; CreateHisto(2);} //8+2
		if(type&32) {varAll.clear();} //reset varAll

		for(int i=0;i<t->GetEntries() ;i++)
			{
			t->GetEntry(i);
			treeVar["ptJet0"]=jetPt[0];
			treeVar["etaJet0"]=jetEta[0];
			treeVar["rhoPF"]=treeVar["rhoIso"];
			treeVar["rhoPFJets"]=treeVar["rho"];
			
			if((treeVar["ptJet0"]<PtMin)||(treeVar["ptJet0"]>PtMax)||(fabs(treeVar["etaJet0"])<EtaMin)||(fabs(treeVar["etaJet0"])>EtaMax)|| (treeVar["rhoPF"]<RhoMin)||(treeVar["rhoPF"]>RhoMax))continue;
			//selection
			if(Invmass < 70 || Invmass > 110)continue;
			if( jetBtag[0] >0.244)continue;
			TLorentzVector mu1,mu2,jet1;
    			mu1.SetPtEtaPhiE(MuPt[0],MuEta[0],MuPhi[0],MuEnergy[0]);
    			mu2.SetPtEtaPhiE(MuPt[1],MuEta[1],MuPhi[1],MuEnergy[1]);
    			double mu12_phi = (mu1+mu2).Phi();
			double muJet_dphi = deltaPhi(jetPhi[0],mu12_phi);
			if(fabs(muJet_dphi)<3.1415-0.5)continue;
			if(jetPt[1]/(mu1+mu2).Pt() >0.3) continue;
				
			//fprintf(stderr,"--C\n");
			
			
			//trigger --only on data
			if( type==1 && !( triggerResult != NULL && triggerResult->size()>1 && triggerResult->at(0) )) continue;
			
			//parton Matching
			double dR_min=999;
			int pos_min=999;
			//int part_min=5;
			
			if(type>1){ //only on MC
			for(int iPart=0;iPart<int(partonPt->size());iPart++)
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
			if(varName=="QGLHisto")treeVar["QGLHisto"] = qgl->computeQGLikelihood2012(jetPt[0],jetEta[0],treeVar["rhoPF"],jetChgPart_QC[0]+jetNeutralPart_ptcut[0]-sub_data,jetPtD_QC[0],jetAxis_QC[1][0]);
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
			{
					int bin=puw->FindBin(treeVar["rho"]);
					float weight=puw->GetBinContent(bin) ;
			treeVar["eventWeight"]=1.;
			treeVar["PUReWeight"]=weight;
			}
		
			//fprintf(stderr,"------G\n");
			if(type&1){
				//printf("passed selection - type 1 --> %.3f - %.3f\n",treeVar[varName],treeVar[varName+"Fwd"]);
				string var=varName;
				alpha=1; beta=0;
					treeVar["eventWeight"]=1.; //data
					treeVar["PUReWeight"]=1;
				FillHisto(h_data,var);
				varAllData.push_back(treeVar[var]);
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
				if( fabs(treeVarInt["pdgIdPartJet0"])== 0) {alpha=a_g;beta=b_g;}
				if( fabs(treeVarInt["pdgIdPartJet0"])== 999) {alpha=a_g;beta=b_g;}

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
					float weight=puw->GetBinContent(bin) ;
					varAll.push_back(TRIO(treeVarInt["pdgIdPartJet0"],treeVar[varName],weight)); //w=-1 default
					//work-around to make it smear all in one shot
					//varAll.push_back(1,treeVar[varName],weight); //w=-1 default
				}
			}
}


void Analyzer::LoadBins(){
		ResetBins();
		PtBins.push_back(  pair<float,float>(20,30) );
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
		EtaBins.push_back(  pair<float,float>(2,2.5) );
		EtaBins.push_back(  pair<float,float>(2.5,3.0) );
		EtaBins.push_back(  pair<float,float>(3,4.7) );
	}



//-----------------------------------DOUBLE MIN FAST ---------------------------------------------------------


//----------------------------------------------------------------------------------------------------



int ComputeDoubleMinZJetHbb(){
	system("[ -f output.root ] && rm output.root");
	TChain *mc=new TChain("Hbb/events");
	TChain *data=new TChain("Hbb/events");
		mc->Add(MCTREE);

		data  ->Add(DATATREE);
	Analyzer A;
	A.nstep=20;
	A.stp0=.01;A.stp1=0.01;
	A.varName="QGLHisto";
	A.CreateHisto();
	A.SetTrees(mc,data);
	A.readStartPoints=1.0;
	fprintf(stderr,"Start Points from DiJet\n");
	//A.filenameStartPoints="../data/SystZJetHbb_2013_07_23.txt" ;
	A.readStartPoints=0;
	//A.filenameStartPoints="../data/SystMB_2013_07_25_DiJetMC.txt";
		freopen("/dev/null","w",stderr);
	fprintf(stderr,"Going to do Span\n");
	A.SpanMin();
	A.varName="QGLMLP";
	//A.SpanMin();
	return 0;
	}

Analyzer *Check(float PtMin,float PtMax,float RhoMin,float RhoMax, float EtaMin,float EtaMax, float alpha, float beta , const char * varName="QGLHisto",float lmin=0,float lmax=1,TCanvas **pC=NULL){
Analyzer *A=new Analyzer();
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add(MCTREE);
data  ->Add(DATATREE);
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
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add(MCTREE);
data  ->Add(DATATREE);
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

void Validation(float PtMin,float PtMax,float RhoMin,float RhoMax, float EtaMin,float EtaMax, const char * varName="QGLHisto",TCanvas **pC=NULL){
Analyzer *A=new Analyzer();
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add(MCTREE);
data  ->Add(DATATREE);
A->varName=varName;
A->RhoMin=RhoMin; A->RhoMax=RhoMax;A->PtMin=PtMin;A->PtMax=PtMax; A->EtaMin=EtaMin;A->EtaMax=EtaMax;
A->nBins=50;A->xMin=0;A->xMax=1.0;
A->CreateHisto();
A->SetTrees(mc,data);

A->Loop(data,1);
A->Loop(mc,32);

TCanvas *c=A->DrawValidation();
if(pC!=NULL)  (*pC)=c;

}

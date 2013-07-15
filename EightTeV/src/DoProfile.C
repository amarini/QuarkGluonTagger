#include "ComputeDoubleMinDiJet.C"

class ProfAnalyzer :public Analyzer{
public:
	ProfAnalyzer():Analyzer(){runMin=0;runMax=500000;};
	bool ExtraCuts(){ 
			//true: event accepted; false event rejected
			if(treeVarInt["runNo"] < 1000 ) return true; //mc
			if (treeVarInt["runNo"]<runMin) return false ; //data
			if (treeVarInt["runNo"]>runMax) return false ;
			return true;};
	float runMin;
	float runMax;
};

ProfAnalyzer *DoProfile(TCanvas **pC=NULL){
ProfAnalyzer *A=new ProfAnalyzer();
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_P6_Dijets.root");
data  ->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012_JetMon_Dijets.root");
A->nstep=20; 
A->lmin=0;
A->lmax=50;
A->nBins=200;
A->varName="multOOB";
A->xMin=0;
A->xMax=50;
A->PtMin=80.;A->PtMax=120.; A->EtaMin=3.0;A->EtaMax=4.7;
A->CreateHisto();
A->SetTrees(mc,data);
A->alpha=1;
A->beta=0;

float RhoBins[]={0,6,8,10,12,14,16,18,20,24,28,32,36,40};
//float RhoBins[]={10,20};
float Mean_data[100];
float Mean_mc[100];
float MeanError_data[100];
float MeanError_mc[100];
float RMS_data[100];
float RMS_mc[100];

int nBins=sizeof(RhoBins)/sizeof(float) - 1;

TCanvas *c0=new TCanvas("c0","c0",800,800);
c0->Divide(3,4);
for(int iBin=0;iBin<nBins;iBin++)
{
	cout<<"Bin "<<iBin<<" of "<<nBins<<endl;
	cout<<"--- "<<RhoBins[iBin]<<"-"<<RhoBins[iBin+1]<<endl;
	A->RhoMin=RhoBins[iBin]; A->RhoMax=RhoBins[iBin+1];
	A->Loop(mc,2);
	A->Loop(data,1);
	Mean_data[iBin]=A->h_data->GetMean();
	Mean_mc[iBin]=A->h_mc->GetMean();
	MeanError_data[iBin]=A->h_data->GetMeanError();
	MeanError_mc[iBin]=A->h_mc->GetMeanError();
	RMS_data[iBin]=A->h_data->GetRMS();
	RMS_mc[iBin]=A->h_mc->GetRMS();
	c0->cd(iBin+1);
	A->h_data->Draw("P");
	A->h_mc->Draw("HIST SAME");

}


gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
TCanvas *c=new TCanvas("c","c",800,800);

TH1F*Hdata=new TH1F("Mult_data","data",nBins,RhoBins);
TH1F*Hmc=new TH1F("Mult_mc","mc",nBins,RhoBins);
TH1F*H_RMS_data=new TH1F("Mult_RMS_data","data",nBins,RhoBins);
TH1F*H_RMS_mc=new TH1F("Mult_RMS_mc","mc",nBins,RhoBins);
for(int iBin=0;iBin<nBins;iBin++)
	{
	Hdata->SetBinContent(iBin+1,Mean_data[iBin]);	
	Hdata->SetBinError(iBin+1,MeanError_data[iBin]);	
	H_RMS_data->SetBinContent(iBin+1,Mean_data[iBin]);	
	H_RMS_data->SetBinError(iBin+1,RMS_data[iBin]);	

	Hmc->SetBinContent(iBin+1,Mean_mc[iBin]);	
	Hmc->SetBinError(iBin+1,MeanError_mc[iBin]);	
	H_RMS_mc->SetBinContent(iBin+1,Mean_mc[iBin]);	
	H_RMS_mc->SetBinError(iBin+1,RMS_mc[iBin]);	
	}


Hdata->SetMarkerStyle(20);
Hdata->SetMarkerColor(kBlack);
Hdata->SetLineColor(kBlack);
H_RMS_data->SetMarkerStyle(0);
H_RMS_data->SetLineColor(kBlack);
H_RMS_data->SetLineStyle(3);
H_RMS_data->SetLineWidth(2);
H_RMS_data->SetFillStyle(0);
H_RMS_data->SetFillColor(kGray+2);

Hmc->SetMarkerStyle(24);
Hmc->SetMarkerColor(54);
Hmc->SetLineColor(54);
H_RMS_mc->SetMarkerStyle(0);
H_RMS_mc->SetLineColor(54);
H_RMS_mc->SetLineStyle(2);
H_RMS_mc->SetLineWidth(2);
H_RMS_mc->SetFillStyle(0);
H_RMS_mc->SetFillColor(54);

H_RMS_data->Draw("P E3");
H_RMS_mc->Draw("P E3 SAME");
Hmc->Draw("P SAME");
Hdata->Draw("P SAME");

	TLegend *L=new TLegend(0.4,.7,.6,.89);
	L->AddEntry(Hdata,"data");
	L->AddEntry(Hmc,"mc");
	L->Draw();
if(pC!=NULL) (*pC)=c;

return A;
}


ProfAnalyzer *DoProfileRun(const char *varName="multOOB",TCanvas **pC=NULL){
ProfAnalyzer *A=new ProfAnalyzer();
TChain *mc=new TChain("Hbb/events");
TChain *data=new TChain("Hbb/events");
mc->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_flatQCD_P6_Dijets.root");
data  ->Add("/afs/cern.ch/work/s/sunil/public/forTom/analysis_data2012_JetMon_Dijets.root");
A->nstep=20; 
A->lmin=0;
A->lmax=50;
A->nBins=200;
A->xMin=0;
A->xMax=50;
A->PtMin=80.;A->PtMax=120.; A->EtaMin=3.0;A->EtaMax=4.7;A->RhoMin=0;A->RhoMax=100;
A->CreateHisto();
A->SetTrees(mc,data);
A->alpha=1;
A->beta=0;

vector<TH1F*> rhoPerRun;

float RunBins[]={ 190000,193000,195000,197000,200000,206000,208000,210000};
//float RhoBins[]={10,20};
float Mean_data[100];
float Mean_mc[100];
float MeanError_data[100];
float MeanError_mc[100];
float RMS_data[100];
float RMS_mc[100];

int nBins=sizeof(RunBins)/sizeof(float) - 1;
//PUWeights
A->varName="rho";
A->nBins=100;
A->xMin=0;
A->xMax=50;
rhoPerRun.reserve(nBins);
TH1F* puw=(TH1F*)A->puw->Clone("default");
A->nBins=puw->GetNbinsX();
A->xMin=puw->GetBinLowEdge(1);
A->xMax=puw->GetBinLowEdge(A->nBins+1);
for(int iBin=0;iBin<nBins;iBin++)
{
	cout<<"Bin "<<iBin<<" of "<<nBins<<endl;
	A->runMin=RunBins[iBin]; A->runMax=RunBins[iBin+1];
	A->Loop(data,1);
	A->Loop(mc,2);
	rhoPerRun.push_back((TH1F*)A->h_data->Clone(Form("puw_%d",iBin)));
	rhoPerRun[iBin]->Divide(A->h_mc);
	rhoPerRun[iBin]->Multiply(puw); //reweighting on top of the default one
	
}
A->varName=varName;
A->nBins=200;
A->xMin=0;
A->xMax=50;

TCanvas *c0=new TCanvas("c0","c0",800,800);
c0->Divide(3,4);
for(int iBin=0;iBin<nBins;iBin++)
{
	cout<<"Bin "<<iBin<<" of "<<nBins<<endl;
	cout<<"--- "<<RunBins[iBin]<<"-"<<RunBins[iBin+1]<<endl;
	A->runMin=RunBins[iBin]; A->runMax=RunBins[iBin+1];
	A->puw=rhoPerRun[iBin];
	A->Loop(mc,2);
	A->Loop(data,1);
	Mean_data[iBin]=A->h_data->GetMean();
	Mean_mc[iBin]=A->h_mc->GetMean();
	MeanError_data[iBin]=A->h_data->GetMeanError();
	MeanError_mc[iBin]=A->h_mc->GetMeanError();
	RMS_data[iBin]=A->h_data->GetRMS();
	RMS_mc[iBin]=A->h_mc->GetRMS();
	c0->cd(iBin+1);
	A->h_data->Draw("P");
	A->h_mc->Draw("HIST SAME");

}


gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
TCanvas *c=new TCanvas("c","c",800,800);

TH1F*Hdata=new TH1F("Mult_data","data",nBins,RunBins);
TH1F*Hmc=new TH1F("Mult_mc","mc",nBins,RunBins);
TH1F*H_RMS_data=new TH1F("Mult_RMS_data","data",nBins,RunBins);
TH1F*H_RMS_mc=new TH1F("Mult_RMS_mc","mc",nBins,RunBins);
for(int iBin=0;iBin<nBins;iBin++)
	{
	Hdata->SetBinContent(iBin+1,Mean_data[iBin]);	
	Hdata->SetBinError(iBin+1,MeanError_data[iBin]);	
	H_RMS_data->SetBinContent(iBin+1,Mean_data[iBin]);	
	H_RMS_data->SetBinError(iBin+1,RMS_data[iBin]);	

	Hmc->SetBinContent(iBin+1,Mean_mc[iBin]);	
	Hmc->SetBinError(iBin+1,MeanError_mc[iBin]);	
	H_RMS_mc->SetBinContent(iBin+1,Mean_mc[iBin]);	
	H_RMS_mc->SetBinError(iBin+1,RMS_mc[iBin]);	
	}


Hdata->SetMarkerStyle(20);
Hdata->SetMarkerColor(kBlack);
Hdata->SetLineColor(kBlack);
H_RMS_data->SetMarkerStyle(0);
H_RMS_data->SetLineColor(kBlack);
H_RMS_data->SetLineStyle(3);
H_RMS_data->SetLineWidth(2);
H_RMS_data->SetFillStyle(0);
H_RMS_data->SetFillColor(kGray+2);

Hmc->SetMarkerStyle(24);
Hmc->SetMarkerColor(54);
Hmc->SetLineColor(54);
H_RMS_mc->SetMarkerStyle(0);
H_RMS_mc->SetLineColor(54);
H_RMS_mc->SetLineStyle(2);
H_RMS_mc->SetLineWidth(2);
H_RMS_mc->SetFillStyle(0);
H_RMS_mc->SetFillColor(54);

H_RMS_data->GetXaxis()->SetTitle("runNum");
H_RMS_data->GetYaxis()->SetTitle("Mult. (mean)");

H_RMS_data->Draw("P E3");
H_RMS_mc->Draw("P E3 SAME");
Hmc->Draw("P SAME");
Hdata->Draw("P SAME");

	TLegend *L=new TLegend(0.4,.7,.6,.89);
	L->AddEntry(Hdata,"data");
	L->AddEntry(Hmc,"mc");
	L->Draw();
if(pC!=NULL) (*pC)=c;

return A;
}


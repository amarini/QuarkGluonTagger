#include <cmath>
#include "TH1D.h"
#include "../interface/QGSyst.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "DrawBase.h"



void drawSinglePlot( DrawBase* db, TH1D* h1_qgl, TH1D* h1_qglSyst, float ptMin, float ptMax, float etaMin, float etaMaxi, float rhoMin, float rhoMax );

int main() {


  TFile* file = TFile::Open("/afs/cern.ch/user/a/amarini/work/GluonTag/ZJet/ZJet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_2.root");
  TTree* tree = (TTree*)file->Get("tree_passedEvents");


  float rho;
  tree->SetBranchAddress("rhoPF", &rho );
  float ptJet;
  tree->SetBranchAddress("ptJet0", &ptJet );
  float etaJet;
  tree->SetBranchAddress("etaJet0", &etaJet );
  float qglJet;
  tree->SetBranchAddress("QGLikelihood2012Jet0", &qglJet );

  QGSyst qgsyst;
  qgsyst.ReadDatabase("../data/SystDatabase.txt");
  qgsyst.SetTagger("QGLHisto");


  TFile* outfile = TFile::Open("prova.root", "RECREATE");
  outfile->cd(); 

  TH1D* h1_qglJet_pt3050_eta02_rho015 = new TH1D("qglJet_pt3050_eta02_rho015", "", 100, 0., 1.0001);
  h1_qglJet_pt3050_eta02_rho015->Sumw2();
  TH1D* h1_qglJetSyst_pt3050_eta02_rho015 = new TH1D("qglJetSyst_pt3050_eta02_rho015", "", 100, 0., 1.0001);
  h1_qglJetSyst_pt3050_eta02_rho015->Sumw2();

  TH1D* h1_qglJet_pt5080_eta02_rho015 = new TH1D("qglJet_pt5080_eta02_rho015", "", 100, 0., 1.0001);
  h1_qglJet_pt5080_eta02_rho015->Sumw2();
  TH1D* h1_qglJetSyst_pt5080_eta02_rho015 = new TH1D("qglJetSyst_pt5080_eta02_rho015", "", 100, 0., 1.0001);
  h1_qglJetSyst_pt5080_eta02_rho015->Sumw2();

  TH1D* h1_qglJet_pt80120_eta02_rho015 = new TH1D("qglJet_pt80120_eta02_rho015", "", 100, 0., 1.0001);
  h1_qglJet_pt80120_eta02_rho015->Sumw2();
  TH1D* h1_qglJetSyst_pt80120_eta02_rho015 = new TH1D("qglJetSyst_pt80120_eta02_rho015", "", 100, 0., 1.0001);
  h1_qglJetSyst_pt80120_eta02_rho015->Sumw2();


  TH1D* h1_qglJet_pt3050_eta35_rho015 = new TH1D("qglJet_pt3050_eta35_rho015", "", 100, 0., 1.0001);
  h1_qglJet_pt3050_eta35_rho015->Sumw2();
  TH1D* h1_qglJetSyst_pt3050_eta35_rho015 = new TH1D("qglJetSyst_pt3050_eta35_rho015", "", 100, 0., 1.0001);
  h1_qglJetSyst_pt3050_eta35_rho015->Sumw2();

  TH1D* h1_qglJet_pt5080_eta35_rho015 = new TH1D("qglJet_pt5080_eta35_rho015", "", 100, 0., 1.0001);
  h1_qglJet_pt5080_eta35_rho015->Sumw2();
  TH1D* h1_qglJetSyst_pt5080_eta35_rho015 = new TH1D("qglJetSyst_pt5080_eta35_rho015", "", 100, 0., 1.0001);
  h1_qglJetSyst_pt5080_eta35_rho015->Sumw2();

  TH1D* h1_qglJet_pt80120_eta35_rho015 = new TH1D("qglJet_pt80120_eta35_rho015", "", 100, 0., 1.0001);
  h1_qglJet_pt80120_eta35_rho015->Sumw2();
  TH1D* h1_qglJetSyst_pt80120_eta35_rho015 = new TH1D("qglJetSyst_pt80120_eta35_rho015", "", 100, 0., 1.0001);
  h1_qglJetSyst_pt80120_eta35_rho015->Sumw2();


  int nentries = tree->GetEntries();

  for( unsigned iEntry = 0; iEntry < nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( rho>15. ) continue;


    if( fabs(etaJet)<2. ) {

      if( ptJet > 30. && ptJet < 50. ) {

        h1_qglJet_pt3050_eta02_rho015->Fill( qglJet );
        h1_qglJetSyst_pt3050_eta02_rho015->Fill( qgsyst.Smear(ptJet, etaJet, rho, qglJet) );

      } else if( ptJet > 50. && ptJet < 80. ) {

        h1_qglJet_pt5080_eta02_rho015->Fill( qglJet );
        h1_qglJetSyst_pt5080_eta02_rho015->Fill( qgsyst.Smear(ptJet, etaJet, rho, qglJet) );

      } else if( ptJet > 80. && ptJet < 120. ) {

        h1_qglJet_pt80120_eta02_rho015->Fill( qglJet );
        h1_qglJetSyst_pt80120_eta02_rho015->Fill( qgsyst.Smear(ptJet, etaJet, rho, qglJet) );

      } 

    } else if( fabs(etaJet)>3. && fabs(etaJet)<5. ) {

      if( ptJet > 30. && ptJet < 50. ) {

        h1_qglJet_pt3050_eta35_rho015->Fill( qglJet );
        h1_qglJetSyst_pt3050_eta35_rho015->Fill( qgsyst.Smear(ptJet, etaJet, rho, qglJet) );

      } else if( ptJet > 50. && ptJet < 80. ) {

        h1_qglJet_pt5080_eta35_rho015->Fill( qglJet );
        h1_qglJetSyst_pt5080_eta35_rho015->Fill( qgsyst.Smear(ptJet, etaJet, rho, qglJet) );

      } else if( ptJet > 80. && ptJet < 120. ) {

        h1_qglJet_pt80120_eta35_rho015->Fill( qglJet );
        h1_qglJetSyst_pt80120_eta35_rho015->Fill( qgsyst.Smear(ptJet, etaJet, rho, qglJet) );

      } 

    } // eta

  }


  DrawBase* db = new DrawBase("provaSyst");
  
  drawSinglePlot( db, h1_qglJet_pt3050_eta02_rho015, h1_qglJetSyst_pt3050_eta02_rho015, 30., 50., 0., 2., 0., 15.);
  drawSinglePlot( db, h1_qglJet_pt5080_eta02_rho015, h1_qglJetSyst_pt5080_eta02_rho015, 50., 80., 0., 2., 0., 15.);
  drawSinglePlot( db, h1_qglJet_pt80120_eta02_rho015, h1_qglJetSyst_pt80120_eta02_rho015, 80., 120., 0., 2., 0., 15.);

  drawSinglePlot( db, h1_qglJet_pt3050_eta35_rho015, h1_qglJetSyst_pt3050_eta35_rho015, 30., 50., 3., 5., 0., 15.);
  drawSinglePlot( db, h1_qglJet_pt5080_eta35_rho015, h1_qglJetSyst_pt5080_eta35_rho015, 50., 80., 3., 5., 0., 15.);
  drawSinglePlot( db, h1_qglJet_pt80120_eta35_rho015, h1_qglJetSyst_pt80120_eta35_rho015, 80., 120., 3., 5., 0., 15.);



  outfile->cd();

  h1_qglJet_pt3050_eta02_rho015->Write();
  h1_qglJetSyst_pt3050_eta02_rho015->Write();
  
  h1_qglJet_pt5080_eta02_rho015->Write();
  h1_qglJetSyst_pt5080_eta02_rho015->Write();
  
  h1_qglJet_pt80120_eta02_rho015->Write();
  h1_qglJetSyst_pt80120_eta02_rho015->Write();

  h1_qglJet_pt3050_eta35_rho015->Write();
  h1_qglJetSyst_pt3050_eta35_rho015->Write();
  
  h1_qglJet_pt5080_eta35_rho015->Write();
  h1_qglJetSyst_pt5080_eta35_rho015->Write();
  
  h1_qglJet_pt80120_eta35_rho015->Write();
  h1_qglJetSyst_pt80120_eta35_rho015->Write();

  outfile->Close();
  


  return 0;
  
}


void drawSinglePlot( DrawBase* db, TH1D* h1_qgl, TH1D* h1_qglSyst, float ptMin, float ptMax, float etaMin, float etaMax, float rhoMin, float rhoMax ) {


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_qgl->Rebin(2);
  h1_qglSyst->Rebin(2);

  h1_qgl->SetLineWidth(2);
  h1_qglSyst->SetLineWidth(2);

  float ymax = h1_qgl->GetMaximum() / h1_qgl->Integral();

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., 1.3*ymax);
  h2_axes->SetXTitle("QG LD");
  h2_axes->SetYTitle("Normalized to Unity");

  h1_qglSyst->SetLineColor(kRed);

  h2_axes->Draw();
  h1_qgl->DrawNormalized("Histo same");
  h1_qglSyst->DrawNormalized("Histo same");

  char legendTitle[500];
  sprintf( legendTitle, "p_{T}(%.0f-%.0f), #eta(%.1f-%.1f), #rho(%.0f-%.0f)", ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax);

  TLegend* legend = new TLegend(0.25, 0.7, 0.75, 0.9, legendTitle);
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_qgl, "Before Smearing", "F" );
  legend->AddEntry( h1_qglSyst, "After Smearing", "F" );
  legend->Draw("same");

  char canvasName[500];
  sprintf( canvasName, "prova_pt%.0f_%.0f_eta%.0f_%.0f_rho%.0f_%.0f.eps", ptMin, ptMax, etaMin, etaMax, rhoMin, rhoMax);

  TPaveText* label_top = db->get_labelTop();
  label_top->Draw("same");

  c1->SaveAs(canvasName);

  delete c1;
  delete legend;
  delete h2_axes;

}

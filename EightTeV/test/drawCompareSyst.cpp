#include "DrawBase.h"





void drawSingleComparison( DrawBase* db, const std::string& selection, bool use_rootChi2 );
void compareToPythia( DrawBase* db, const std::string& selection, const std::string& smearingName );



int main() {


  DrawBase* db = new  DrawBase("compareSyst");

  drawSingleComparison( db, "ZJet", true );
  drawSingleComparison( db, "DiJet", true );

  drawSingleComparison( db, "ZJet", false );
  drawSingleComparison( db, "DiJet", false );


  compareToPythia( db, "DiJet", "noQsmear" );
  compareToPythia( db, "DiJet", "nominal" );
  compareToPythia( db, "DiJet", "smearQLikeP6" );
  compareToPythia( db, "DiJet", "noQsmearPt80" );


  return 0;

}



void drawSingleComparison( DrawBase* db, const std::string& selection, bool use_rootChi2 ) {

  TFile* file_nominal = TFile::Open(Form("checkSyst_doubleMin_%s_qglJet_Hpp.root", selection.c_str()));
  TFile* file_noQSmear = TFile::Open(Form("checkSyst_doubleMin_noQsmear_%s_qglJet_Hpp.root", selection.c_str()));
  TFile* file_smearQLikeP6 = TFile::Open(Form("checkSyst_doubleMin_smearQLikeP6_%s_qglJet_Hpp.root", selection.c_str()));
  TFile* file_newTest = TFile::Open(Form("checkSyst_doubleMin_newTest_%s_qglJet_Hpp.root", selection.c_str()));

  std::string suffix = (use_rootChi2) ? "root" : "";

  TGraph* gr_nominal      = (TGraph*)file_nominal      ->Get(Form("chi2%s_vs_pt_eta02", suffix.c_str()));
  TGraph* gr_noQSmear     = (TGraph*)file_noQSmear     ->Get(Form("chi2%s_vs_pt_eta02", suffix.c_str()));
  TGraph* gr_smearQLikeP6 = (TGraph*)file_smearQLikeP6 ->Get(Form("chi2%s_vs_pt_eta02", suffix.c_str()));
  TGraph* gr_newTest      = (TGraph*)file_newTest      ->Get(Form("chi2%s_vs_pt_eta02", suffix.c_str()));

  gr_nominal      ->SetMarkerSize(1.6);
  gr_noQSmear     ->SetMarkerSize(1.6);
  gr_smearQLikeP6 ->SetMarkerSize(1.6);
  gr_newTest      ->SetMarkerSize(1.6);

  gr_nominal      ->SetMarkerStyle(20);
  gr_noQSmear     ->SetMarkerStyle(21);
  gr_smearQLikeP6 ->SetMarkerStyle(24);
  gr_newTest      ->SetMarkerStyle(25);

  gr_nominal      ->SetMarkerColor(46);
  gr_noQSmear     ->SetMarkerColor(38);
  gr_smearQLikeP6 ->SetMarkerColor(29);
  gr_newTest      ->SetMarkerColor(39);


  float yMax =  (use_rootChi2) ? 20. : 40.;

  TH2D* h2_axes = new TH2D( "axes", "", 10, 26., 200., 10, 0., yMax );
  h2_axes->SetXTitle( "Jet p_{T} [GeV]" );
  h2_axes->SetYTitle( "#chi^{2} / NDF" );

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();
  h2_axes->Draw();


  TPaveText* selectionText = new TPaveText( 0.2, 0.8, 0.4, 0.9, "brNDC" );
  selectionText->SetTextSize( 0.038 );
  selectionText->SetFillColor( 0 );
  selectionText->AddText( selection.c_str() );
  selectionText->Draw("same");

  TLegend* legend = new TLegend( 0.5, 0.7, 0.9, 0.9 );
  legend->SetTextSize( 0.038 );
  legend->SetFillColor(kWhite);
  legend->AddEntry( gr_nominal, "PAS", "P" );
  legend->AddEntry( gr_noQSmear, "No quark smearing", "P" );
  legend->AddEntry( gr_smearQLikeP6, "P6 quark smearing", "P" );
  legend->AddEntry( gr_newTest, "New test", "P" );
  legend->Draw("same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  gr_noQSmear     ->Draw("P same");
  gr_smearQLikeP6 ->Draw("P same");
  gr_newTest      ->Draw("P same");
  gr_nominal      ->Draw("P same");

  gPad->RedrawAxis();

  std::string canvasName = "systCompare_chi2";
  if( use_rootChi2 ) canvasName += "root";
  canvasName = canvasName + "_" + selection;

  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";

  c1->SaveAs( canvasName_eps.c_str() );
  c1->SaveAs( canvasName_png.c_str() );

  delete h2_axes;
  delete c1;
  delete legend;

}



void compareToPythia( DrawBase* db, const std::string& selection, const std::string& smearingName ) {


  TFile* file_pythia = TFile::Open(Form("checkSyst_doubleMin_%s_qglJet.root", selection.c_str()));
  TFile* file_herwig;
  if( smearingName != "nominal" ) 
    file_herwig= TFile::Open(Form("checkSyst_doubleMin_%s_%s_qglJet_Hpp.root", smearingName.c_str(), selection.c_str()));
  else
    file_herwig= TFile::Open(Form("checkSyst_doubleMin_%s_qglJet_Hpp.root", selection.c_str()));


  TH1D* h1_eff_quark_pythia = (TH1D*)file_pythia->Get("effNum_centr_quark_syst_thresh1");
  TH1D* h1_eff_quark_herwig = (TH1D*)file_herwig->Get("effNum_centr_quark_syst_thresh1");

  TH1D* h1_eff_gluon_pythia = (TH1D*)file_pythia->Get("effNum_centr_gluon_syst_thresh1");
  TH1D* h1_eff_gluon_herwig = (TH1D*)file_herwig->Get("effNum_centr_gluon_syst_thresh1");

  h1_eff_quark_pythia->SetMarkerSize(1.6);
  h1_eff_gluon_pythia->SetMarkerSize(1.6);
  h1_eff_quark_herwig->SetMarkerSize(1.6);
  h1_eff_gluon_herwig->SetMarkerSize(1.6);

  h1_eff_quark_pythia->SetMarkerStyle(20);
  h1_eff_gluon_pythia->SetMarkerStyle(21);
  h1_eff_quark_herwig->SetMarkerStyle(24);
  h1_eff_gluon_herwig->SetMarkerStyle(25);

  h1_eff_quark_pythia->SetMarkerColor(38);
  h1_eff_gluon_pythia->SetMarkerColor(46);
  h1_eff_quark_herwig->SetMarkerColor(38);
  h1_eff_gluon_herwig->SetMarkerColor(46);


  TCanvas* c1 = new TCanvas("c1", "", 600, 600 );
  c1->cd();
  c1->SetLogx();


  TH2D* h2_axes = new TH2D("axes", "", 10, 30., 520., 10, 0., 1.);
  h2_axes->SetXTitle("Jet p_{T} [GeV]");
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->SetYTitle("Efficiency (LD > 0.5)");
  h2_axes->Draw();

  TPaveText* labelTop = db->get_labelTop();
  TPaveText* labelEta = new TPaveText(0.22, 0.2, 0.36, 0.25, "brNDC");
  labelEta->AddText("|#eta| < 2");
  labelEta->SetTextSize(0.038);
  labelEta->SetFillColor(0);

  TPaveText* labelSelection = new TPaveText(.55,.2,.85,.25,"brNDC");
  labelSelection->AddText(Form("Smeared %s MC",selection.c_str()));
  labelSelection->SetTextSize(0.038);
  labelSelection->SetFillColor(0);


  TLegend* legend_quark = new TLegend( 0.17, 0.79, 0.45, 0.91 );
  legend_quark->SetTextSize(0.04);
  legend_quark->SetFillColor(kWhite);
  legend_quark->AddEntry( h1_eff_quark_pythia, "Quark (Pythia6)", "P");
  legend_quark->AddEntry( h1_eff_quark_herwig, "Quark (Herwig++)", "P");

  TLegend* legend_gluon = new TLegend( 0.55, 0.79, 0.83, 0.91 );
  legend_gluon->SetTextSize(0.04);
  legend_gluon->SetFillColor(kWhite);
  legend_gluon->AddEntry( h1_eff_gluon_pythia, "Gluon (Pythia6)", "P");
  legend_gluon->AddEntry( h1_eff_gluon_herwig, "Gluon (Herwig++)", "P");

  legend_gluon->Draw("same");
  legend_quark->Draw("same");
  labelEta->Draw("same");
  labelTop->Draw("same");
  labelSelection->Draw("same");
  gStyle->SetErrorX(0);

  h1_eff_quark_pythia->Draw("P same");
  h1_eff_gluon_pythia->Draw("P same");
  h1_eff_quark_herwig->Draw("P same");
  h1_eff_gluon_herwig->Draw("P same");

  gPad->RedrawAxis();

  std::string canvasName = "compareToPythia_" + smearingName + "_" + selection;
  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";

  c1->SaveAs( canvasName_eps.c_str());
  c1->SaveAs( canvasName_png.c_str());

  delete c1;
  delete legend_quark;
  delete legend_gluon;
  delete h2_axes;

}


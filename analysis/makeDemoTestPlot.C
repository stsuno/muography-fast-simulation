//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Understanding primitve variables.
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "utils.h"

void alabel(int itype) {
  if (itype == 1) {
    TPad *tp = new TPad("","",0.5,0.85,0.95,0.95,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.5);
    text -> DrawLatex(0.25,0.5,"(A)");
    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(1); line->SetLineStyle(1); line->DrawLine(0.05,0.55,0.2,0.55);
  }
  if (itype == 2) {
    TPad *tp = new TPad("","",0.65,0.8,0.92,0.92,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.34);
    text -> DrawLatex(0.2,0.60,"w/ anomaly");
    text -> DrawLatex(0.2,0.10,"w/o anomaly");
    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(2); line->SetLineStyle(1); line->DrawLine(0.0,0.7,0.16,0.7);
    line->SetLineColor(1); line->SetLineStyle(1); line->DrawLine(0.0,0.2,0.16,0.2);
  }
}

void makeDemoTestPlot(){

  char xtit[50];
  char ytit[50];
  const int nfil = 2;
  std::string filename[nfil];

  std::string plottit = "plot.ps";
  filename[0] = "hist_test_anomaly.root";
  filename[1] = "hist_test_noanomaly.root";
 
  TH1F* hist_track_nhits[nfil];
  TH1F* hist_track_chi2[nfil];
  TH1F* hist_track_C[nfil];
  TH1F* hist_track_D[nfil];
  TH1F* hist_track_errA[nfil];
  TH1F* hist_track_errB[nfil];
  TH1F* hist_track_errC[nfil];
  TH1F* hist_track_errD[nfil];
  TH1F* hist_track_phi[nfil];
  TH1F* hist_track_theta[nfil];
  TH1F* hist_muon_energy[nfil];
  TH1F* hist_muon_loss[nfil];
  TH1F* hist_muon_detenergy[nfil];
  TH1F* hist_muon_charge[nfil];
  TH1F* hist_muon_phi[nfil];
  TH1F* hist_muon_theta[nfil];
  TH1F* hist_diff_phi[nfil];
  TH1F* hist_diff_theta[nfil];
  TH1F* hist_residual_x1[nfil];
  TH1F* hist_residual_x2[nfil];
  TH1F* hist_residual_x3[nfil];
  TH1F* hist_residual_x4[nfil];
  TH1F* hist_residual_y1[nfil];
  TH1F* hist_residual_y2[nfil];
  TH1F* hist_residual_y3[nfil];
  TH1F* hist_residual_y4[nfil];
  TH2F* hist2D_xy50[nfil];
  TH2F* hist2D_xy100[nfil];

  TFile *ffl[nfil];
  for (int i = 0; i < nfil; i++) {
    ffl[i] = new TFile(filename[i].c_str(),"read");

    hist_track_nhits[i] = (TH1F*) ffl[i]->Get("hist_track_nhits");
    hist_track_chi2[i] = (TH1F*) ffl[i]->Get("hist_track_chi2");
    hist_track_C[i] = (TH1F*) ffl[i]->Get("hist_track_C");
    hist_track_D[i] = (TH1F*) ffl[i]->Get("hist_track_D");
    hist_track_errA[i] = (TH1F*) ffl[i]->Get("hist_track_errA");
    hist_track_errB[i] = (TH1F*) ffl[i]->Get("hist_track_errB");
    hist_track_errC[i] = (TH1F*) ffl[i]->Get("hist_track_errC");
    hist_track_errD[i] = (TH1F*) ffl[i]->Get("hist_track_errD");
    hist_track_phi[i] = (TH1F*) ffl[i]->Get("hist_track_phi");
    hist_track_theta[i] = (TH1F*) ffl[i]->Get("hist_track_theta");
    hist_muon_energy[i] = (TH1F*) ffl[i]->Get("hist_muon_energy");
    hist_muon_loss[i] = (TH1F*) ffl[i]->Get("hist_muon_loss");
    hist_muon_detenergy[i] = (TH1F*) ffl[i]->Get("hist_muon_detenergy");
    hist_muon_charge[i] = (TH1F*) ffl[i]->Get("hist_muon_charge");
    hist_muon_phi[i] = (TH1F*) ffl[i]->Get("hist_muon_phi");
    hist_muon_theta[i] = (TH1F*) ffl[i]->Get("hist_muon_theta");
    hist_diff_phi[i] = (TH1F*) ffl[i]->Get("hist_diff_phi");
    hist_diff_theta[i] = (TH1F*) ffl[i]->Get("hist_diff_theta");
    hist_residual_x1[i] = (TH1F*) ffl[i]->Get("hist_residual_x1");
    hist_residual_x2[i] = (TH1F*) ffl[i]->Get("hist_residual_x2");
    hist_residual_x3[i] = (TH1F*) ffl[i]->Get("hist_residual_x3");
    hist_residual_x4[i] = (TH1F*) ffl[i]->Get("hist_residual_x4");
    hist_residual_y1[i] = (TH1F*) ffl[i]->Get("hist_residual_y1");
    hist_residual_y2[i] = (TH1F*) ffl[i]->Get("hist_residual_y2");
    hist_residual_y3[i] = (TH1F*) ffl[i]->Get("hist_residual_y3");
    hist_residual_y4[i] = (TH1F*) ffl[i]->Get("hist_residual_y4");
    hist2D_xy50[i] = (TH2F*) ffl[i]->Get("hist2D_xy50");
    hist2D_xy100[i] = (TH2F*) ffl[i]->Get("hist2D_xy100");

  }

  TPostScript *ps = new TPostScript(plottit.c_str(),-111);

  TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);
  c1->SetBorderSize(0);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBorderMode(0);
  gStyle -> SetPadBorderSize(0);
  gStyle -> SetCanvasBorderMode(0);

  gStyle -> SetPadLeftMargin(0.15);
  //*gStyle -> SetPadRightMargin(0.05);
  //*gStyle -> SetPadTopMargin(0.05);
  gStyle -> SetPadBottomMargin(0.15);

  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.5);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelOffset(0.003,"X");
  gStyle->SetLabelOffset(0.003,"Y");

  gStyle->SetPalette(1);
  gStyle->SetNdivisions(510);

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

//  gStyle -> SetPadLeftMargin(0.17);
  gStyle -> SetPadLeftMargin(0.18);
//  gStyle->SetTitleYOffset(1.5);
  gStyle->SetTitleYOffset(2.0);
//  gStyle -> SetPadRightMargin(0.114);

  ps -> NewPage();

  //========
  // page.1
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_nhits[0]->GetSum(),
             hist_track_nhits[1]->GetSum(),
             hist_track_nhits[0],
             hist_track_nhits[1],
             0,0,"Number of hits","Events");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_chi2[0]->GetSum(),
             hist_track_chi2[1]->GetSum(),
             hist_track_chi2[0],
             hist_track_chi2[1],
             0,0,"Track #chi^{2}/ndf","Events");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_C[0]->GetSum(),
             hist_track_C[1]->GetSum(),
             hist_track_C[0],
             hist_track_C[1],
             0,0,"Detector X[cm]","Events");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_D[0]->GetSum(),
             hist_track_D[1]->GetSum(),
             hist_track_D[0],
             hist_track_D[1],
             0,0,"Detector Y[cm]","Events");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.2
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_errA[0]->GetSum(),
             hist_track_errA[1]->GetSum(),
             hist_track_errA[0],
             hist_track_errA[1],
             0,0,"Track Error A","Events");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_errB[0]->GetSum(),
             hist_track_errB[1]->GetSum(),
             hist_track_errB[0],
             hist_track_errB[1],
             0,0,"Track Error B","Events");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_errC[0]->GetSum(),
             hist_track_errC[1]->GetSum(),
             hist_track_errC[0],
             hist_track_errC[1],
             0,0,"Track Error C","Events");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_errD[0]->GetSum(),
             hist_track_errD[1]->GetSum(),
             hist_track_errD[0],
             hist_track_errD[1],
             0,0,"Track Error D","Events");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.3
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_phi[0]->GetSum(),
             hist_track_phi[1]->GetSum(),
             hist_track_phi[0],
             hist_track_phi[1],
             0,0,"Track #phi","Events");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_track_theta[0]->GetSum(),
             hist_track_theta[1]->GetSum(),
             hist_track_theta[0],
             hist_track_theta[1],
             0,0,"Track #theta","Events");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_diff_phi[0]->GetSum(),
             hist_diff_phi[1]->GetSum(),
             hist_diff_phi[0],
             hist_diff_phi[1],
             0,0,"#Delta#phi","Events");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_diff_theta[0]->GetSum(),
             hist_diff_theta[1]->GetSum(),
             hist_diff_theta[0],
             hist_diff_theta[1],
             0,0,"#Delta#theta","Events");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.4
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_muon_energy[0]->GetSum(),
             hist_muon_energy[1]->GetSum(),
             hist_muon_energy[0],
             hist_muon_energy[1],
             0,0,"Muon energy[GeV]","Events");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_muon_loss[0]->GetSum(),
             hist_muon_loss[1]->GetSum(),
             hist_muon_loss[0],
             hist_muon_loss[1],
             0,0,"Muon energy loss[MeV]","Events");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_muon_detenergy[0]->GetSum(),
             hist_muon_detenergy[1]->GetSum(),
             hist_muon_detenergy[0],
             hist_muon_detenergy[1],
             0,0,"Muon energy at detector[GeV]","Events");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_muon_charge[0]->GetSum(),
             hist_muon_charge[1]->GetSum(),
             hist_muon_charge[0],
             hist_muon_charge[1],
             0,0,"Muon charge","Events");
//  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.5
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_x1[0]->GetSum(),
             hist_residual_x1[1]->GetSum(),
             hist_residual_x1[0],
             hist_residual_x1[1],
             0,0,"Residual L1X","Events");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_x2[0]->GetSum(),
             hist_residual_x2[1]->GetSum(),
             hist_residual_x2[0],
             hist_residual_x2[1],
             0,0,"Residual L2X","Events");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_x3[0]->GetSum(),
             hist_residual_x3[1]->GetSum(),
             hist_residual_x3[0],
             hist_residual_x3[1],
             0,0,"Residual L3X","Events");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_x4[0]->GetSum(),
             hist_residual_x4[1]->GetSum(),
             hist_residual_x4[0],
             hist_residual_x4[1],
             0,0,"Residual L4X","Events");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.6
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_y1[0]->GetSum(),
             hist_residual_y1[1]->GetSum(),
             hist_residual_y1[0],
             hist_residual_y1[1],
             0,0,"Residual L1Y","Events");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_y2[0]->GetSum(),
             hist_residual_y2[1]->GetSum(),
             hist_residual_y2[0],
             hist_residual_y2[1],
             0,0,"Residual L2Y","Events");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_y3[0]->GetSum(),
             hist_residual_y3[1]->GetSum(),
             hist_residual_y3[0],
             hist_residual_y3[1],
             0,0,"Residual L3Y","Events");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,2,1,1,1,
             hist_residual_y4[0]->GetSum(),
             hist_residual_y4[1]->GetSum(),
             hist_residual_y4[0],
             hist_residual_y4[1],
             0,0,"Residual L4Y","Events");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.7
  //========

  int nbinX = hist2D_xy50[0]->GetNbinsX();
  int nbinY = hist2D_xy50[0]->GetNbinsY();
  double xmin = hist2D_xy50[0]->GetXaxis()->GetXmin();
  double xmax = hist2D_xy50[0]->GetXaxis()->GetXmax();
  double ymin = hist2D_xy50[0]->GetYaxis()->GetXmin();
  double ymax = hist2D_xy50[0]->GetYaxis()->GetXmax();
  TH2F* hist2D_sub50 = new TH2F("hist2D_sub50","",nbinX,xmin,xmax,nbinY,ymin,ymax);
  TH2F* hist2D_sig50 = new TH2F("hist2D_sig50","",nbinX,xmin,xmax,nbinY,ymin,ymax);

  TH2F* hist2D_sub100 = new TH2F("hist2D_sub100","",nbinX,xmin,xmax,nbinY,ymin,ymax);
  TH2F* hist2D_sig100 = new TH2F("hist2D_sig100","",nbinX,xmin,xmax,nbinY,ymin,ymax);

  TH2F* hist2D_sub200 = new TH2F("hist2D_sub200","",nbinX,xmin,xmax,nbinY,ymin,ymax);
  TH2F* hist2D_sig200 = new TH2F("hist2D_sig200","",nbinX,xmin,xmax,nbinY,ymin,ymax);

  for (int i=0; i<nbinX; i++) {
    for (int j=0; j<nbinY; j++) {
      double ch1 = hist2D_xy50[0]->GetBinContent(i+1,j+1);
      double ch2 = hist2D_xy50[1]->GetBinContent(i+1,j+1);
      hist2D_sub50 -> SetBinContent(i+1,j+1,ch1-ch2);
      hist2D_sig50 -> SetBinContent(i+1,j+1,(ch1-ch2)*(ch1-ch2)/(ch1+ch2));
      double ci1 = hist2D_xy100[0]->GetBinContent(i+1,j+1);
      double ci2 = hist2D_xy100[1]->GetBinContent(i+1,j+1);
      hist2D_sub100 -> SetBinContent(i+1,j+1,ci1-ci2);
      hist2D_sig100 -> SetBinContent(i+1,j+1,(ci1-ci2)*(ci1-ci2)/(ci1+ci2));
//      std::cout << i+1 << " " << j+1 << " " << (ch1-ch2)*(ch1-ch2)/(ch1+ch2) << std::endl;
    }
  }


  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=50)[cm]", "y(z=50)[cm]");
  hist2D_xy50[0] -> Draw("colz,same");
  c1->cd(2);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=50)[cm]", "y(z=50)[cm]");
  hist2D_xy50[1] -> Draw("colz,same");
  c1->cd(3);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=50)[cm]", "y(z=50)[cm]");
//  hist2D_sub50 -> SetMaximum(2.0);
//  hist2D_sub50 -> SetMinimum(-2.0);
  hist2D_sub50 -> Draw("colz,same");
  c1->cd(4);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=50)[cm]", "y(z=50)[cm]");
  hist2D_sig50 -> SetMaximum(2.0);
  hist2D_sig50 -> SetMinimum(0.0);
  hist2D_sig50 -> Draw("colz,same");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.8
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=100)[cm]", "y(z=100)[cm]");
  hist2D_xy100[0] -> Draw("colz,same");
  c1->cd(2);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=100)[cm]", "y(z=100)[cm]");
  hist2D_xy100[1] -> Draw("colz,same");
  c1->cd(3);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=100)[cm]", "y(z=100)[cm]");
//  hist2D_sub100 -> SetMaximum(2.0);
//  hist2D_sub100 -> SetMinimum(-2.0);
  hist2D_sub100 -> Draw("colz,same");
  c1->cd(4);   c1->SetLogy(0);
  waku(-50.0, 50.0, -50.0, 50.0, "x(z=100)[cm]", "y(z=100)[cm]");
  hist2D_sig100 -> SetMaximum(2.0);
  hist2D_sig100 -> SetMinimum(0.0);
  hist2D_sig100 -> Draw("colz,same");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.9
  //========

  for (int i = 0; i < nfil; i++) {

    hist_track_nhits[i] -> Clear();
    hist_track_chi2[i] -> Clear();
    hist_track_C[i] -> Clear();
    hist_track_D[i] -> Clear();
    hist_track_errA[i] -> Clear();
    hist_track_errB[i] -> Clear();
    hist_track_errC[i] -> Clear();
    hist_track_errD[i] -> Clear();
    hist_track_phi[i] -> Clear();
    hist_track_theta[i] -> Clear();
    hist_muon_energy[i] -> Clear();
    hist_muon_loss[i] -> Clear();
    hist_muon_detenergy[i] -> Clear();
    hist_muon_charge[i] -> Clear();
    hist_muon_phi[i] -> Clear();
    hist_muon_theta[i] -> Clear();
    hist_diff_phi[i] -> Clear();
    hist_diff_theta[i] -> Clear();
    hist_residual_x1[i] -> Clear();
    hist_residual_x2[i] -> Clear();
    hist_residual_x3[i] -> Clear();
    hist_residual_x4[i] -> Clear();
    hist_residual_y1[i] -> Clear();
    hist_residual_y2[i] -> Clear();
    hist_residual_y3[i] -> Clear();
    hist_residual_y4[i] -> Clear();
    hist2D_xy50[i] -> Clear();
    hist2D_xy100[i] -> Clear();


    ffl[i] -> Close();
  }

  ps -> Close();
}


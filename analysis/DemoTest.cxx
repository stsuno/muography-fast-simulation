#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "DemoTest.h"
//#include "Helix.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <math.h>

//class Point3D;
//class Helix;

void DemoTest::initialize(TString rout) {
  m_rootFile = new TFile(rout,"recreate");

  hist_track_nhits = new TH1F("hist_track_nhits","", 10,   0.5,  10.5); 
  hist_track_chi2  = new TH1F("hist_track_chi2", "", 20,   0.0,   0.4); 
  hist_track_C     = new TH1F("hist_track_C",    "", 50, -50.0,  50.0); 
  hist_track_D     = new TH1F("hist_track_D",    "", 50, -50.0,  50.0); 
  hist_track_errA  = new TH1F("hist_track_errA", "", 20,   0.0,   0.2); 
  hist_track_errB  = new TH1F("hist_track_errB", "", 20,   0.0,   0.2); 
  hist_track_errC  = new TH1F("hist_track_errC", "", 50,   0.0,   2.0); 
  hist_track_errD  = new TH1F("hist_track_errD", "", 50,   0.0,   2.0); 
  hist_track_phi   = new TH1F("hist_track_phi",  "",100,-180.0, 180.0); 
  hist_track_theta = new TH1F("hist_track_theta","",100,   0.0,  90.0); 

  hist_muon_energy    = new TH1F("hist_muon_energy",   "",100, 0.0, 10.0);
  hist_muon_loss      = new TH1F("hist_muon_loss",     "",100, 0.0, 2000.0);
  hist_muon_detenergy = new TH1F("hist_muon_detenergy","",100, 0.0, 10.0);
  hist_muon_charge    = new TH1F("hist_muon_charge",   "",  2,-1.5, 1.5);
  hist_muon_phi       = new TH1F("hist_muon_phi",      "", 90,-180.0,180.0);
  hist_muon_theta     = new TH1F("hist_muon_theta",    "", 90, 0.0,  90.0);
  hist_diff_phi       = new TH1F("hist_diff_phi",      "",100,-5.0, 5.0);
  hist_diff_theta     = new TH1F("hist_diff_theta",    "",100,-5.0, 5.0);

  hist_residual_x1 = new TH1F("hist_residual_x1","",100,-5.0, 5.0);
  hist_residual_x2 = new TH1F("hist_residual_x2","",100,-5.0, 5.0);
  hist_residual_x3 = new TH1F("hist_residual_x3","",100,-5.0, 5.0);
  hist_residual_x4 = new TH1F("hist_residual_x4","",100,-5.0, 5.0);
  hist_residual_y1 = new TH1F("hist_residual_y1","",100,-5.0, 5.0);
  hist_residual_y2 = new TH1F("hist_residual_y2","",100,-5.0, 5.0);
  hist_residual_y3 = new TH1F("hist_residual_y3","",100,-5.0, 5.0);
  hist_residual_y4 = new TH1F("hist_residual_y4","",100,-5.0, 5.0);

  hist2D_xy50  = new TH2F("hist2D_xy50", "",50, -50, 50,50, -50, 50);
  hist2D_xy100 = new TH2F("hist2D_xy100","",50, -50, 50,50, -50, 50);

}

void DemoTest::execute() {

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << "Total Number of Evetns = " << nentries << std::endl;

  Int_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (jentry%5000==0) { std::cout << jentry << " events processed..." << std::endl; }

    //================
    // Fill histogram
    //================
    hist_track_nhits -> Fill(1.0*num_hits);
    hist_track_chi2  -> Fill(reduced_chisq);
    hist_track_C     -> Fill(C);
    hist_track_D     -> Fill(D);
    hist_track_errA  -> Fill(errorA);
    hist_track_errB  -> Fill(errorB);
    hist_track_errC  -> Fill(errorC);
    hist_track_errD  -> Fill(errorD);
    hist_track_phi   -> Fill(phi);
    hist_track_theta -> Fill(theta);

    hist_muon_energy    -> Fill(muon_energy/1000.0);
    if (muon_loss>0.0) {
      hist_muon_loss      -> Fill(muon_loss);
    }
    hist_muon_detenergy -> Fill((muon_energy-muon_loss)/1000.0);
    hist_muon_charge    -> Fill(1.0*muon_charge);
    hist_muon_phi       -> Fill(muon_phi);
    hist_muon_theta     -> Fill(muon_theta);
    hist_diff_phi       -> Fill(phi-muon_phi);
    hist_diff_theta     -> Fill(theta-muon_theta);

    hist_residual_x1 -> Fill((*residual_x)[0]*10.0);
    hist_residual_x2 -> Fill((*residual_x)[1]*10.0);
    hist_residual_x3 -> Fill((*residual_x)[2]*10.0);
    hist_residual_x4 -> Fill((*residual_x)[3]*10.0);
    hist_residual_y1 -> Fill((*residual_y)[0]*10.0);
    hist_residual_y2 -> Fill((*residual_y)[1]*10.0);
    hist_residual_y3 -> Fill((*residual_y)[2]*10.0);
    hist_residual_y4 -> Fill((*residual_y)[3]*10.0);

//    if (theta<30.0) {
    hist2D_xy50  -> Fill(A*50.0+C, B*50.0+D);
    hist2D_xy100 -> Fill(A*100.0+C,B*100.0+D);
//    }

  }
}

void DemoTest::finalize() {
  m_rootFile->Write();
  m_rootFile->Close();
}



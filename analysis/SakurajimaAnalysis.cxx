#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "SakurajimaAnalysis.h"
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

void SakurajimaAnalysis::initialize(TString rout) {
  m_rootFile = new TFile(rout,"recreate");

  hist_track_nhits = new TH1F("hist_track_nhits","", 10,   5.5,  15.5); 
  hist_track_chi2  = new TH1F("hist_track_chi2", "", 100,   0.0,   5.0); 
//  hist_track_C     = new TH1F("hist_track_C",    "", 50, -50.0,  50.0); 
//  hist_track_D     = new TH1F("hist_track_D",    "", 50, -50.0,  50.0); 
  hist_track_C     = new TH1F("hist_track_C",    "", 50, -100.0,  100.0); 
  hist_track_D     = new TH1F("hist_track_D",    "", 50, -100.0,  100.0); 
  hist_track_phi   = new TH1F("hist_track_phi",  "",100,-180.0, 180.0); 
  hist_track_theta = new TH1F("hist_track_theta","",100,   0.0,  90.0); 

  hist2D_mountain  = new TH2F("hist2D_mountain","",80,-200.0,200.0,40,200.0,400.0); 
  hist2D_sky       = new TH2F("hist2D_sky",     "",80,-200.0,200.0,40,200.0,400.0); 

}

void SakurajimaAnalysis::execute() {

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
    hist_track_chi2  -> Fill(chisq/(2.0*num_hits-4.0));
    hist_track_C     -> Fill(C);
    hist_track_D     -> Fill(D);
    hist_track_phi   -> Fill(phi);
    hist_track_theta -> Fill(theta);

    // position at mountain
    double distanceZ = 2.5*1000.0*100.0;   // [cm]
    double offsetY = 200.0;                // [m]
    double wirespacing = 1.25;

    if (phi<0.0) {
      hist2D_mountain -> Fill(wirespacing*(-A*distanceZ+C)/1000.0, wirespacing*(-B*distanceZ+D)/1000.0+offsetY);
    }
    else {
      hist2D_sky -> Fill(wirespacing*(A*distanceZ+C)/1000.0, wirespacing*(B*distanceZ+D)/1000.0+offsetY);
    }

  }
}

void SakurajimaAnalysis::finalize() {
  m_rootFile->Write();
  m_rootFile->Close();
}



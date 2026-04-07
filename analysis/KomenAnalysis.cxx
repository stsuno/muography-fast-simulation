#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "KomenAnalysis.h"
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

void KomenAnalysis::initialize(TString rout) {
  m_rootFile = new TFile(rout,"recreate");

  hist_track_theta = new TH1F("hist_track_theta","",100,   0.0,  90.0); 
  hist_track_phi   = new TH1F("hist_track_phi",  "",100,-180.0, 180.0); 
}

void KomenAnalysis::execute() {

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
    hist_track_theta -> Fill(theta);
    hist_track_phi   -> Fill(phi);

  }
}

void KomenAnalysis::finalize() {
  m_rootFile->Write();
  m_rootFile->Close();
}



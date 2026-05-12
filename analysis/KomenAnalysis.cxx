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
  hist_narrow_theta = new TH1F("hist_narrow_theta","",20, 0.0, 10.0); 

  hist_mb_phi    = new TH1F("hist_mb_phi",  "", 45,-180.0, 180.0); 
  hist_mb_theta  = new TH1F("hist_mb_theta","", 32,   0.0,   8.0); 
  diff_mb_phi    = new TH1F("diff_mb_phi",   "",100, -5.0, 5.0);
  diff_mb_theta  = new TH1F("diff_mb_theta", "", 33, -4.5, 4.5);
  diff_mb_deltaR = new TH1F("diff_mb_deltaR","",100,  0.0,10.0);
  diff_mb_C      = new TH1F("diff_mb_C",     "",100,-200.0,200.0);
  diff_mb_D      = new TH1F("diff_mb_D",     "",100,-200.0,200.0);
  diff_mb_dist   = new TH1F("diff_mb_dist",  "", 25,  0.0,250.0);

  hist_mc_phi    = new TH1F("hist_mc_phi",  "", 45,-180.0, 180.0); 
  hist_mc_theta  = new TH1F("hist_mc_theta","", 24,  22.0,  28.0); 
  diff_mc_phi    = new TH1F("diff_mc_phi",   "",100, -5.0, 5.0);
  diff_mc_theta  = new TH1F("diff_mc_theta", "", 33, -4.5, 4.5);
  diff_mc_deltaR = new TH1F("diff_mc_deltaR","",100,  0.0,10.0);
  diff_mc_C      = new TH1F("diff_mc_C",     "",100,-200.0,200.0);
  diff_mc_D      = new TH1F("diff_mc_D",     "",100,-200.0,200.0);
  diff_mc_dist   = new TH1F("diff_mc_dist",  "", 25,  0.0,250.0);

  hist_ma_phi    = new TH1F("hist_ma_phi",  "",  40, -8.0,  8.0); 
  hist_ma_theta  = new TH1F("hist_ma_theta","",  48, 37.0, 43.0); 
  diff_ma_phi    = new TH1F("diff_ma_phi",   "", 40, -4.0, 4.0);
  diff_ma_theta  = new TH1F("diff_ma_theta", "", 21, -2.5, 2.5);
  diff_ma_deltaR = new TH1F("diff_ma_deltaR","", 50,  0.0, 4.0);
  diff_ma_C      = new TH1F("diff_ma_C",     "",100,-200.0,200.0);
  diff_ma_D      = new TH1F("diff_ma_D",     "",100,-100.0,100.0);
  diff_ma_dist   = new TH1F("diff_ma_dist",  "", 25,  0.0,250.0);

  hist_bc_phi    = new TH1F("hist_bc_phi",  "",  40, -0.05, 0.05); 
  hist_bc_theta  = new TH1F("hist_bc_theta","",  12, 74.0, 80.0); 

  hist_mb_simpletheta = new TH1F("hist_mb_simpletheta","", 9, -0.5,  8.5); 
  hist_mc_simpletheta = new TH1F("hist_mc_simpletheta","", 7, 21.5, 28.5); 
  hist_ma_simpletheta = new TH1F("hist_ma_simpletheta","", 7, 36.5, 43.5); 
  hist_bc_simpletheta = new TH1F("hist_bc_simpletheta","", 7, 73.5, 80.5); 

}

void KomenAnalysis::execute() {

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << "Total Number of Evetns = " << nentries << std::endl;

  //======================================================
  TChain* tree2 = new TChain("track_tree");
  track_tree *chktrk;
  Long64_t ntot = 0;
  int matched_event = -1;
  bool flag_matched = true;
  if (m_mode==2) {
    TString rin = m_FDirectoryName + m_checkfile + "/track_tree";
    std::cout << "Checking File name  : " << m_checkfile << std::endl;
    tree2->AddFile(rin,-1); 

    ntot = tree2->GetEntriesFast();
    std::cout << "Total Number of events in the check file : " << ntot << std::endl;

    TTree *tnt = tree2;
    chktrk = new track_tree(tnt);
  }
  //======================================================

  Int_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (jentry%5000==0) { std::cout << jentry << " events processed..." << std::endl; }

    float theta_chk = -10000.0;
    float phi_chk   = -10000.0;
    float A_chk     = -10000.0;
    float B_chk     = -10000.0;
    float C_chk     = -10000.0;
    float D_chk     = -10000.0;

    if (m_mode==2) {
      int nb2 = 0;
      flag_matched = false;
      int indx_matched = -1;
      for (int i=matched_event+1; i<ntot; i++) {
        Long64_t ichk = tree2->LoadTree(i);
        if (ichk<0) { break; }
        nb2 = tree2->GetEntry(i);
        if (chktrk->event_number==event_number) {
          flag_matched = true;
          indx_matched = i;
          break;
        }
      }
      if (indx_matched>-1) {
        matched_event = indx_matched;
      }
    }

    if (flag_matched) {
      theta_chk = (chktrk->theta);
      phi_chk   = (chktrk->phi);
      A_chk     = (chktrk->A);
      B_chk     = (chktrk->B);
      C_chk     = (chktrk->C);
      D_chk     = (chktrk->D);

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
      hist_narrow_theta -> Fill(theta);

      hist_mb_phi    -> Fill(0.5*(phi+phi_chk));
      hist_mb_theta  -> Fill(0.5*(theta+theta_chk));
      diff_mb_phi    -> Fill(phi-phi_chk);
      diff_mb_theta  -> Fill(theta-theta_chk);
      diff_mb_deltaR -> Fill(std::sqrt((phi-phi_chk)*(phi-phi_chk)+(theta-theta_chk)*(theta-theta_chk)));

      double X1 = C_chk+m_positionMWPC.X();
      double Y1 = D_chk+m_positionMWPC.Y();
      double X2 = A*m_positionMWPC.Z()+C;
      double Y2 = B*m_positionMWPC.Z()+D;

      diff_mb_C      -> Fill(X1-X2);
      diff_mb_D      -> Fill(Y1-Y2);
      diff_mb_dist   -> Fill(std::sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)));
      hist_mb_simpletheta -> Fill(0.5*(theta+theta_chk));

      hist_mc_phi    -> Fill(0.5*(phi+phi_chk));
      hist_mc_theta  -> Fill(0.5*(theta+theta_chk));
      diff_mc_phi    -> Fill(phi-phi_chk);
      diff_mc_theta  -> Fill(theta-theta_chk);
      diff_mc_deltaR -> Fill(std::sqrt((phi-phi_chk)*(phi-phi_chk)+(theta-theta_chk)*(theta-theta_chk)));
      diff_mc_C      -> Fill(X1-X2);
      diff_mc_D      -> Fill(Y1-Y2);
      diff_mc_dist   -> Fill(std::sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)));
      hist_mc_simpletheta -> Fill(0.5*(theta+theta_chk));

      double phi_ma = 90.0-phi;
      double theta_ma = 90.0-theta;

      double expX_ma = -A_chk*m_positionMWPC.Z()+C_chk;
      double expY_ma = -D_chk;

      hist_ma_phi    -> Fill(0.5*(phi_ma+phi_chk));
      hist_ma_theta  -> Fill(0.5*(theta_ma+theta_chk));
      diff_ma_phi    -> Fill(phi_ma-phi_chk);
      diff_ma_theta  -> Fill(theta_ma-theta_chk);
      diff_ma_deltaR -> Fill(std::sqrt((phi_ma-phi_chk)*(phi_ma-phi_chk)+(theta_ma-theta_chk)*(theta_ma-theta_chk)));
      diff_ma_C      -> Fill(C-expX_ma);
      diff_ma_D      -> Fill(D-expY_ma);
      diff_ma_dist   -> Fill(std::sqrt((C-expX_ma)*(C-expX_ma)+(D-expY_ma)*(D-expY_ma)));
      hist_ma_simpletheta -> Fill(0.5*(theta_ma+theta_chk));

      hist_bc_phi    -> Fill(phi);
      hist_bc_theta  -> Fill(theta);
      hist_bc_simpletheta -> Fill(theta);
    }

  }
}

void KomenAnalysis::finalize() {
  m_rootFile->Write();
  m_rootFile->Close();
}



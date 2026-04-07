#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "makeTree.h"
//#include "Helix.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <math.h>

using namespace std;

//class Point3D;
//class Helix;

void makeTree::initialize(TString rout) {
  m_rootFile = new TFile(rout,"recreate");

  m_trktree = new TTree("trktree","trktree");
  m_trktree->SetDirectory(m_rootFile);
  m_trktree->Branch("runNumber",&m_runNumber);
  m_trktree->Branch("eventNumber",&m_eventNumber);
  m_trktree->Branch("lumiBlock",&m_lumiB);
  m_trktree->Branch("averagePU",&m_averagePU);
  m_trktree->Branch("eventPU",&m_eventPU);
  m_trktree->Branch("mcFlag",&m_mcFlag);
  m_trktree->Branch("pileupWeight",&m_pileupWeight);
  m_trktree->Branch("mcWeight",&m_mcWeight);

  m_trktree->Branch("pVtxX",&m_pVtxX);
  m_trktree->Branch("pVtxY",&m_pVtxY);
  m_trktree->Branch("pVtxZ",&m_pVtxZ);

  m_trktree->Branch("pVtxXErr",&m_pVtxXErr);
  m_trktree->Branch("pVtxYErr",&m_pVtxYErr);
  m_trktree->Branch("pVtxZErr",&m_pVtxZErr);

  m_trktree->Branch("truthPVtxX",&m_truthPVtxX);
  m_trktree->Branch("truthPVtxY",&m_truthPVtxY);
  m_trktree->Branch("truthPVtxZ",&m_truthPVtxZ);

  m_trktree->Branch("trackNPixelHits",&m_nPixelHits);
  m_trktree->Branch("trackNPixelHoles",&m_nPixelHoles);
  m_trktree->Branch("trackNPixelLayers",&m_nPixelLayers);
  m_trktree->Branch("trackNPixelOutliers",&m_nPixelOutliers);
  m_trktree->Branch("trackNIBLHits",&m_nIBLHits);
  m_trktree->Branch("trackNSCTHits",&m_nSCTHits);
  m_trktree->Branch("trackNSCTHoles",&m_nSCTHoles);
  m_trktree->Branch("trackNTRTHits",&m_nTRTHits);
  m_trktree->Branch("trackNBLHits",&m_nBLHits);
  m_trktree->Branch("trackNSplitIBLHits",&m_nIBLSplitHits);
  m_trktree->Branch("trackNSplitBLHits",&m_nBLSplitHits);
  m_trktree->Branch("trackNSharedIBLHits",&m_nIBLSharedHits);
  m_trktree->Branch("trackNSharedBLHits",&m_nBLSharedHits);
  m_trktree->Branch("trackNL2Hits",&m_nL2Hits);
  m_trktree->Branch("trackNL3Hits",&m_nL3Hits);
  m_trktree->Branch("trackNECHits",&m_nECHits);

  m_trktree->Branch("trackNIBLExpectedHits",&m_nIBLExpectedHits);
  m_trktree->Branch("trackNBLExpectedHits",&m_nBLExpectedHits);
  m_trktree->Branch("trackNSCTSharedHits",&m_nSCTSharedHits);

  m_trktree->Branch("trackCharge",&m_trackCharge);
  m_trktree->Branch("trackPt",&m_trackPt);
  m_trktree->Branch("trackPhi",&m_trackPhi);
  m_trktree->Branch("trackEta",&m_trackEta);
  m_trktree->Branch("trackTheta",&m_trackTheta);
  m_trktree->Branch("trackqOverP",&m_trackqOverP);
  m_trktree->Branch("trackD0",&m_trackD0);
  m_trktree->Branch("trackZ0",&m_trackZ0);
  m_trktree->Branch("trackD0Err",&m_trackD0Err);
  m_trktree->Branch("trackZ0Err",&m_trackZ0Err);
  m_trktree->Branch("trackqOverPErr",&m_trackqOverPErr);
  m_trktree->Branch("trackDeltaZSinTheta",&m_trackDeltaZSinTheta);

  m_trktree->Branch("trackClass",&m_trackClass);
  m_trktree->Branch("trackPassCut",&m_trackPassCut);

  m_trktree->Branch("trackPixeldEdx",&m_pixeldEdx);
  m_trktree->Branch("trackNPixeldEdx",&m_nPixeldEdx);

  m_trktree->Branch("trueD0",&m_trueD0);
  m_trktree->Branch("trueZ0",&m_trueZ0);
  m_trktree->Branch("truePhi",&m_truePhi);
  m_trktree->Branch("trueTheta",&m_trueTheta);
  m_trktree->Branch("trueqOverP",&m_trueqOverP);
  m_trktree->Branch("truePdgId",&m_truepdgid);
  m_trktree->Branch("trueBarcode",&m_truebarcode);
  m_trktree->Branch("truthMatchProb",&m_truthMatchProb);

  m_trktree->Branch("genVtxR",&m_genVtxR);
  m_trktree->Branch("genVtxZ",&m_genVtxZ);
  m_trktree->Branch("parentFlavour",&m_parentFlav);

  m_trktree->Branch("minTrkdR",&m_mindR);
  m_trktree->Branch("hitSplitIBL",&m_hitSplitIBL);
  m_trktree->Branch("hitSplitBL",&m_hitSplitBL);
  m_trktree->Branch("hitSplitL2",&m_hitSplitL2);
  m_trktree->Branch("hitSplitL3",&m_hitSplitL3);
 
  m_trktree->Branch("hitIsEndCap",&m_hitIsEndCap);
  m_trktree->Branch("hitIsHole",&m_hitIsHole);
  m_trktree->Branch("hitIsOutlier",&m_hitIsOutlier);
  m_trktree->Branch("hitIs3D",&m_hitIs3D);
  m_trktree->Branch("hitIsIBL",&m_hitIsIBL);
  m_trktree->Branch("hitLayer",&m_hitLayer);
  m_trktree->Branch("hitCharge",&m_hitCharge);
  m_trktree->Branch("hitToT",&m_hitToT);
  m_trktree->Branch("hitLVL1A",&m_hitLVL1A);
  m_trktree->Branch("hitNPixelX",&m_hitNPixelX);
  m_trktree->Branch("hitNPixelY",&m_hitNPixelY);
  m_trktree->Branch("hitEtaModule",&m_hitEtaModule);
  m_trktree->Branch("hitPhiModule",&m_hitPhiModule);
  m_trktree->Branch("hitVBias",&m_hitVBias);
  m_trktree->Branch("hitVDep",&m_hitVDep);
  m_trktree->Branch("hitTemp",&m_hitTemp);
  m_trktree->Branch("hitLorentzShift",&m_hitLorentzShift);
  m_trktree->Branch("hitIsSplit",&m_hitIsSplit);
  m_trktree->Branch("hitGlobalX",&m_hitGlobalX);
  m_trktree->Branch("hitGlobalY",&m_hitGlobalY);
  m_trktree->Branch("hitGlobalZ",&m_hitGlobalZ);
  m_trktree->Branch("hitLocalX",&m_hitLocalX);
  m_trktree->Branch("hitLocalY",&m_hitLocalY);

  m_trktree->Branch("g4LocalX",&m_g4LocalX);
  m_trktree->Branch("g4LocalY",&m_g4LocalY);
  m_trktree->Branch("g4Barcode",&m_g4Barcode);
  m_trktree->Branch("g4PdgId",&m_g4PdgId);
  m_trktree->Branch("g4EnergyDeposit",&m_g4EnergyDeposit);

  m_trktree->Branch("trkLocalX",&m_trkLocalX);
  m_trktree->Branch("trkLocalY",&m_trkLocalY);
  m_trktree->Branch("hitLocalErrorX",&m_hitLocalErrorX);
  m_trktree->Branch("hitLocalErrorY",&m_hitLocalErrorY);
  m_trktree->Branch("trkLocalErrorX",&m_trkLocalErrorX);
  m_trktree->Branch("trkLocalErrorY",&m_trkLocalErrorY);
  m_trktree->Branch("unbiasedResidualX",&m_unbiasedResidualX);
  m_trktree->Branch("unbiasedResidualY",&m_unbiasedResidualY);
  m_trktree->Branch("biasedResidualX",&m_biasedResidualX);
  m_trktree->Branch("biasedResidualY",&m_biasedResidualY);
  m_trktree->Branch("unbiasedPullX",&m_unbiasedPullX);
  m_trktree->Branch("unbiasedPullY",&m_unbiasedPullY);
  m_trktree->Branch("biasedPullX",&m_biasedPullX);
  m_trktree->Branch("biasedPullY",&m_biasedPullY);

  m_trktree->Branch("trk2G4DistX",&m_minTrkG4DistX);
  m_trktree->Branch("trk2G4DistY",&m_minTrkG4DistY);

  m_trktree->Branch("trk2DistX",&m_minTrkDistX);
  m_trktree->Branch("trk2DistY",&m_minTrkDistY);
  m_trktree->Branch("trkPhiOnSurface",&m_trkPhiOnSurface);
  m_trktree->Branch("trkThetaOnSurface",&m_trkThetaOnSurface);

  m_trktree->Branch("hitNContributingPtcs",&m_clus_nContributingPtcs);
  m_trktree->Branch("hitNContributingPU",&m_clus_nContributingPU);
  m_trktree->Branch("hitContributingPtcsBarcode",&m_clus_contributingPtcsBarcode);

  m_trktree->Branch("pixelIBLIndex",&m_clus_pixelIBLIndex);
  m_trktree->Branch("pixelIBLRow",&m_clus_pixelIBLRow);
  m_trktree->Branch("pixelIBLCol",&m_clus_pixelIBLCol);
  m_trktree->Branch("pixelIBLCharge",&m_clus_pixelIBLCharge);
  m_trktree->Branch("pixelIBLToT",&m_clus_pixelIBLToT);

  m_trktree->Branch("pixelBLIndex",&m_clus_pixelBLIndex);
  m_trktree->Branch("pixelBLRow",&m_clus_pixelBLRow);
  m_trktree->Branch("pixelBLCol",&m_clus_pixelBLCol);
  m_trktree->Branch("pixelBLCharge",&m_clus_pixelBLCharge);
  m_trktree->Branch("pixelBLToT",&m_clus_pixelBLToT);

  m_trktree->Branch("pixelL2Index",&m_clus_pixelL2Index);
  m_trktree->Branch("pixelL2Row",&m_clus_pixelL2Row);
  m_trktree->Branch("pixelL2Col",&m_clus_pixelL2Col);
  m_trktree->Branch("pixelL2Charge",&m_clus_pixelL2Charge);
  m_trktree->Branch("pixelL2ToT",&m_clus_pixelL2ToT);

  m_trktree->Branch("pixelL3Index",&m_clus_pixelL3Index);
  m_trktree->Branch("pixelL3Row",&m_clus_pixelL3Row);
  m_trktree->Branch("pixelL3Col",&m_clus_pixelL3Col);
  m_trktree->Branch("pixelL3Charge",&m_clus_pixelL3Charge);
  m_trktree->Branch("pixelL3ToT",&m_clus_pixelL3ToT);

}

void makeTree::execute() {

  Long64_t nentries = fChain->GetEntriesFast();

  cout << "Total Number of Evetns = " << nentries << endl;

  Int_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (jentry%1000==0) { cout << jentry << " events processed..." << endl; }

    //================
    // Fill histogram
    //================
    m_runNumber = runNumber;

    m_eventNumber = eventNumber;
    m_lumiB       = lumiBlock;
    m_averagePU   = averagePU;
    m_eventPU     = eventPU;
    m_mcFlag      = mcFlag;
    m_pileupWeight = pileupWeight;
    m_mcWeight     = mcWeight;

    m_pVtxX = pVtxX;
    m_pVtxY = pVtxY;
    m_pVtxZ = pVtxZ;

    m_pVtxXErr = pVtxXErr;
    m_pVtxYErr = pVtxYErr;
    m_pVtxZErr = pVtxZErr;

    m_truthPVtxX = truthPVtxX;
    m_truthPVtxY = truthPVtxY;
    m_truthPVtxZ = truthPVtxZ;

    int numTracks = (int)(*trackPt).size();
    for (int i=0; i<numTracks; i++) {
      if ((*trackPt)[i]<2.0) { continue; }
      if ((*trackDeltaZSinTheta)[i]>3.0) { continue; }

      m_trackPt             = (*trackPt)[i];
      m_trackPhi            = (*trackPhi)[i];
      m_trackEta            = (*trackEta)[i];
      m_trackTheta          = (*trackTheta)[i];
      m_trackCharge         = (*trackCharge)[i];
      m_trackqOverP         = (*trackqOverP)[i];
      m_trackD0             = (*trackD0)[i];
      m_trackZ0             = (*trackZ0)[i];
      m_trackD0Err          = (*trackD0Err)[i];
      m_trackZ0Err          = (*trackZ0Err)[i];
      m_trackqOverPErr      = (*trackqOverPErr)[i];
      m_trackDeltaZSinTheta = (*trackDeltaZSinTheta)[i];

      m_trackClass          = (*trackClass)[i];
      m_trackPassCut        = (*trackPassCut)[i];

      m_nPixelHits          = (*trackNPixelHits)[i];
      m_nPixelHoles         = (*trackNPixelHoles)[i];
      m_nPixelLayers        = (*trackNPixelLayers)[i];
      m_nPixelOutliers      = (*trackNPixelOutliers)[i];
      m_nIBLHits            = (*trackNIBLUsedHits)[i];
      m_nSCTHits            = (*trackNSCTHits)[i];
      m_nSCTHoles           = (*trackNSCTHoles)[i];
      m_nTRTHits            = (*trackNTRTHits)[i];
      m_nBLHits             = (*trackNBLUsedHits)[i];
      m_nIBLSplitHits       = (*trackNSplitIBLHits)[i];
      m_nBLSplitHits        = (*trackNSplitBLHits)[i];
      m_nIBLSharedHits      = (*trackNSharedIBLHits)[i];
      m_nBLSharedHits       = (*trackNSharedBLHits)[i];
      m_nL2Hits             = (*trackNL1Hits)[i];
      m_nL3Hits             = (*trackNL2Hits)[i];
      m_nECHits             = (*trackNECHits)[i];

      m_nIBLExpectedHits    = (*trackNIBLExpectedHits)[i];
      m_nBLExpectedHits     = (*trackNBLExpectedHits)[i];
      m_nSCTSharedHits      = (*trackNSCTSharedHits)[i];

      m_pixeldEdx           = (*trackTruncateddEdx)[i];
      m_nPixeldEdx          = (*trackNTruncateddEdx)[i];

      m_trueD0              = (*trued0)[i];
      m_trueZ0              = (*truez0)[i];
      m_truePhi             = (*truephi0)[i];
      m_trueTheta           = (*truetheta)[i];
      m_trueqOverP          = (*trueqoverp)[i];
      m_truepdgid           = (*truepdgid)[i];
      m_truebarcode         = (*truebarcode)[i];
      m_truthMatchProb      = (*truthmatchprob)[i];

      m_genVtxR             = (*genVtxR)[i];
      m_genVtxZ             = (*genVtxZ)[i];
      m_parentFlav          = (*parentFlavour)[i];

      m_mindR               = (*trackMindR)[i];

      m_hitIsEndCap.clear();
      m_hitIsHole.clear();
      m_hitIsOutlier.clear();
      m_hitIs3D.clear();
      m_hitIsIBL.clear();
      m_hitLayer.clear();
      m_hitCharge.clear();
      m_hitToT.clear();
      m_hitLVL1A.clear();
      m_hitNPixelX.clear();
      m_hitNPixelY.clear();
      m_hitEtaModule.clear();
      m_hitPhiModule.clear();
      m_hitVBias.clear();
      m_hitVDep.clear();
      m_hitTemp.clear();
      m_hitLorentzShift.clear();
      m_hitIsSplit.clear();
      m_hitGlobalX.clear();
      m_hitGlobalY.clear();
      m_hitGlobalZ.clear();
      m_hitLocalX.clear();
      m_hitLocalY.clear();
      m_g4LocalX.clear();
      m_g4LocalY.clear();
      m_g4Barcode.clear();
      m_g4PdgId.clear();
      m_g4EnergyDeposit.clear();
      m_trkLocalX.clear();
      m_trkLocalY.clear();
      m_hitLocalErrorX.clear();
      m_hitLocalErrorY.clear();
      m_trkLocalErrorX.clear();
      m_trkLocalErrorY.clear();
      m_unbiasedResidualX.clear();
      m_unbiasedResidualY.clear();
      m_biasedResidualX.clear();
      m_biasedResidualY.clear();
      m_unbiasedPullX.clear();
      m_unbiasedPullY.clear();
      m_biasedPullX.clear();
      m_biasedPullY.clear();
      m_minTrkG4DistX.clear();
      m_minTrkG4DistY.clear();
      m_minTrkDistX.clear();
      m_minTrkDistY.clear();
      m_trkPhiOnSurface.clear();
      m_trkThetaOnSurface.clear();
      m_clus_nContributingPtcs.clear();
      m_clus_nContributingPU.clear();
      m_clus_contributingPtcsBarcode.clear();
      m_clus_pixelIBLIndex.clear();
      m_clus_pixelIBLRow.clear();
      m_clus_pixelIBLCol.clear();
      m_clus_pixelIBLCharge.clear();
      m_clus_pixelIBLToT.clear();
      m_clus_pixelBLIndex.clear();
      m_clus_pixelBLRow.clear();
      m_clus_pixelBLCol.clear();
      m_clus_pixelBLCharge.clear();
      m_clus_pixelBLToT.clear();
      m_clus_pixelL2Index.clear();
      m_clus_pixelL2Row.clear();
      m_clus_pixelL2Col.clear();
      m_clus_pixelL2Charge.clear();
      m_clus_pixelL2ToT.clear();
      m_clus_pixelL3Index.clear();
      m_clus_pixelL3Row.clear();
      m_clus_pixelL3Col.clear();
      m_clus_pixelL3Charge.clear();
      m_clus_pixelL3ToT.clear();

      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {

        int isIBL = 0;
        int is3D = 0;
        if ((*hitIsEndCap)[i][j]==0 && (*hitLayer)[i][j]==0) {  // IBL
          isIBL = 1;
          m_hitSplitIBL = (*hitIsSplit)[i][j];
          if ((*hitEtaModule)[i][j]>5)  { is3D=1; }
          if ((*hitEtaModule)[i][j]<-6) { is3D=1; }

          if ((*RdoCharge)[i][j].size()>0) {
            for (int k=0; k<(int)(*RdoCharge)[i][j].size(); k++) {
              m_clus_pixelIBLIndex.push_back(k);
              m_clus_pixelIBLCharge.push_back((*RdoCharge)[i][j][k]);
              m_clus_pixelIBLToT.push_back((*RdoToT)[i][j][k]);
              m_clus_pixelIBLRow.push_back((*RdoPhi)[i][j][k]);
              m_clus_pixelIBLCol.push_back((*RdoEta)[i][j][k]);
            }
          }
        }
        else if ((*hitIsEndCap)[i][j]==0 && (*hitLayer)[i][j]==1) {  // B-layer
          m_hitSplitBL  = (*hitIsSplit)[i][j];

          if ((*RdoCharge)[i][j].size()>0) {
            for (int k=0; k<(int)(*RdoCharge)[i][j].size(); k++) {
              m_clus_pixelBLIndex.push_back(k);
              m_clus_pixelBLCharge.push_back((*RdoCharge)[i][j][k]);
              m_clus_pixelBLToT.push_back((*RdoToT)[i][j][k]);
              m_clus_pixelBLRow.push_back((*RdoPhi)[i][j][k]);
              m_clus_pixelBLCol.push_back((*RdoEta)[i][j][k]);
            }
          }
        }
        else if ((*hitIsEndCap)[i][j]==0 && (*hitLayer)[i][j]==2) {  // Layer-1
          m_hitSplitL2  = (*hitIsSplit)[i][j];

          if ((*RdoCharge)[i][j].size()>0) {
            for (int k=0; k<(int)(*RdoCharge)[i][j].size(); k++) {
              m_clus_pixelL2Index.push_back(k);
              m_clus_pixelL2Charge.push_back((*RdoCharge)[i][j][k]);
              m_clus_pixelL2ToT.push_back((*RdoToT)[i][j][k]);
              m_clus_pixelL2Row.push_back((*RdoPhi)[i][j][k]);
              m_clus_pixelL2Col.push_back((*RdoEta)[i][j][k]);
            }
          }
        }
        else if ((*hitIsEndCap)[i][j]==0 && (*hitLayer)[i][j]==3) {  // Layer-2
          m_hitSplitL3  = (*hitIsSplit)[i][j];

          if ((*RdoCharge)[i][j].size()>0) {
            for (int k=0; k<(int)(*RdoCharge)[i][j].size(); k++) {
              m_clus_pixelL3Index.push_back(k);
              m_clus_pixelL3Charge.push_back((*RdoCharge)[i][j][k]);
              m_clus_pixelL3ToT.push_back((*RdoToT)[i][j][k]);
              m_clus_pixelL3Row.push_back((*RdoPhi)[i][j][k]);
              m_clus_pixelL3Col.push_back((*RdoEta)[i][j][k]);
            }
          }
        }


        m_hitIsEndCap.push_back((*hitIsEndCap)[i][j]);

        if ((*hitIsHole)[i].size()==0) {
          m_hitIsHole.push_back(0);
        }
        else {
          for (int k=0; k<(int)(*hitIsHole)[i].size(); k++) {
            if ((*hitIsHole)[i][k]==(*hitLayer)[i][j]) {
              m_hitIsHole.push_back(1);
              break;
            }
          }
        }
//        m_hitIsOutlier.push_back((*hitIsOutlier)[i][j]);
        m_hitIsOutlier.push_back(0);

        m_hitIs3D.push_back(is3D);
        m_hitIsIBL.push_back(isIBL);
        m_hitLayer.push_back((*hitLayer)[i][j]);
        m_hitCharge.push_back((*hitCharge)[i][j]);
        m_hitToT.push_back((*hitToT)[i][j]);
        m_hitLVL1A.push_back((*hitLVL1A)[i][j]);
        m_hitNPixelX.push_back((*hitNPixelX)[i][j]);
        m_hitNPixelY.push_back((*hitNPixelY)[i][j]);
        m_hitEtaModule.push_back((*hitEtaModule)[i][j]);
        m_hitPhiModule.push_back((*hitPhiModule)[i][j]);
        m_hitVBias.push_back((*hitVBias)[i][j]);
        m_hitVDep.push_back(0.0);
        m_hitTemp.push_back((*hitTemp)[i][j]);
        m_hitLorentzShift.push_back((*hitLorentzShift)[i][j]);
        m_hitIsSplit.push_back((*hitIsSplit)[i][j]);
        m_hitGlobalX.push_back((*hitGlobalX)[i][j]);
        m_hitGlobalY.push_back((*hitGlobalY)[i][j]);
        m_hitGlobalZ.push_back((*hitGlobalZ)[i][j]);
        m_hitLocalX.push_back((*hitLocalX)[i][j]);
        m_hitLocalY.push_back((*hitLocalY)[i][j]);
        m_trkLocalX.push_back((*trkLocalX)[i][j]);
        m_trkLocalY.push_back((*trkLocalY)[i][j]);
        m_hitLocalErrorX.push_back(0.0);
        m_hitLocalErrorY.push_back(0.0);
        m_trkLocalErrorX.push_back(0.0);
        m_trkLocalErrorY.push_back(0.0);
        m_unbiasedResidualX.push_back((*unbiasedResidualX)[i][j]);
        m_unbiasedResidualY.push_back((*unbiasedResidualY)[i][j]);
        m_biasedResidualX.push_back(0.0);
        m_biasedResidualY.push_back(0.0);
        m_unbiasedPullX.push_back((*unbiasedPullX)[i][j]);
        m_unbiasedPullY.push_back((*unbiasedPullY)[i][j]);
        m_biasedPullX.push_back(0.0);
        m_biasedPullY.push_back(0.0);
        m_trkPhiOnSurface.push_back((*trkPhiOnSurface)[i][j]);
        m_trkThetaOnSurface.push_back((*trkThetaOnSurface)[i][j]);

        // Fill G4 information
        if (mcFlag && (*truebarcode)[i]!=-1 && (*siHitBarcode)[i].size()) {
          int indexSiHit = -1;

          for (int k=0; k<(int)(*siHitBarcode)[i][j].size(); k++) {
            if ((*siHitBarcode)[i][j][k]==(*truebarcode)[i]) { indexSiHit=k; break; } 
          }
          if (indexSiHit>-1) {
            m_g4LocalX.push_back(((*siHitEndPosX)[i][j][indexSiHit]+(*siHitStartPosX)[i][j][indexSiHit])/2.0);
            m_g4LocalY.push_back(((*siHitEndPosY)[i][j][indexSiHit]+(*siHitStartPosY)[i][j][indexSiHit])/2.0);
            m_g4Barcode.push_back((*siHitBarcode)[i][j][indexSiHit]);
            m_g4PdgId.push_back((*siHitPdgId)[i][j][indexSiHit]);
            m_g4EnergyDeposit.push_back((*siHitEnergyDeposit)[i][j][indexSiHit]);
          }
          else {
            m_g4LocalX.push_back(-999.0);
            m_g4LocalY.push_back(-999.0);
            m_g4Barcode.push_back(-999.0);
            m_g4PdgId.push_back(-999.0);
            m_g4EnergyDeposit.push_back(-999.0);
          }

          if (indexSiHit>-1) {
            std::vector<int> otherParticles;  
            for (int k=0; k<(int)(*siHitBarcode)[i][j].size(); k++) {
              if (k==indexSiHit) { continue; }
              if ((*siHitBarcode)[i][j][k]>0 && (*siHitBarcode)[i][j][k]<200000) { otherParticles.push_back(k); }
            }
            if (otherParticles.size()>0) {
              float minDeltaR = 100000.0;
              float minDeltaX = 999.0;
              float minDeltaY = 999.0;

              for (int k=0; k<(int)otherParticles.size(); k++) {
                float refX = ((*siHitEndPosX)[i][j][indexSiHit]+(*siHitStartPosX)[i][j][indexSiHit])/2.0;
                float refY = ((*siHitEndPosY)[i][j][indexSiHit]+(*siHitStartPosY)[i][j][indexSiHit])/2.0;
                float localX = ((*siHitEndPosX)[i][j][otherParticles[k]]+(*siHitStartPosX)[i][j][otherParticles[k]])/2.0;
                float localY = ((*siHitEndPosY)[i][j][otherParticles[k]]+(*siHitStartPosY)[i][j][otherParticles[k]])/2.0;
                float dx = refX-localX;
                float dy = refY-localY;
                float deltaR = TMath::Sqrt(dx*dx+dy*dy);
                if (deltaR<minDeltaR) {
                  minDeltaR = deltaR;
                  minDeltaX = TMath::Abs(dx);
                  minDeltaY = TMath::Abs(dy);
                }
              }
              m_minTrkG4DistX.push_back(minDeltaX);
              m_minTrkG4DistY.push_back(minDeltaY);
            }
            else {
              m_minTrkG4DistX.push_back(999.0);
              m_minTrkG4DistY.push_back(999.0);
            }
          }
        }

        // look for neighbor tracks
        float minDeltaTrk = 100000.0;
        float minDeltaX = 999.0;
        float minDeltaY = 999.0;
        for (int k=0; k<numTracks; k++) {
          if (k==i) { continue; }

          float dEta = (*trackEta)[k]-(*trackEta)[i];
          float dPhi = (*trackPhi)[k]-(*trackPhi)[i];
          float dr = TMath::Sqrt(dEta*dEta+dPhi*dPhi);

          if ((*trackPt)[k]>0.5 && dr<0.85) {
            for (int l=0; l<(int)(*hitLayer)[k].size(); l++) {
              int thisBec   = (*hitIsEndCap)[k][l];
              int thisLayer = (*hitLayer)[k][l];
              if (thisBec==0) {
                if (thisLayer==(*hitIsEndCap)[i][j] && (*hitEtaModule)[k][l]==(*hitEtaModule)[i][j] && (*hitPhiModule)[k][l]==(*hitPhiModule)[i][j]) {
                  float thisX = (*hitLocalX)[k][l]-(*unbiasedResidualX)[k][l];
                  float thisY = (*hitLocalY)[k][l]-(*unbiasedResidualY)[k][l];
                  float refX  = (*hitLocalX)[i][j]-(*unbiasedResidualX)[i][j];
                  float refY  = (*hitLocalY)[i][j]-(*unbiasedResidualY)[i][j];
                  float distTrkX = TMath::Abs(thisX-refX);
                  float distTrkY = TMath::Abs(thisY-refY);
                  if (TMath::Sqrt(distTrkX*distTrkX+distTrkY*distTrkY)<minDeltaTrk && TMath::Sqrt(distTrkX*distTrkX+distTrkY*distTrkY)>2e-3) {
                    minDeltaTrk = sqrt(distTrkX*distTrkX+distTrkY*distTrkY);
                    minDeltaX = distTrkX;
                    minDeltaY = distTrkY;
                  }
                }
              }
            }
          }
        }
        m_minTrkDistX.push_back(minDeltaX);		
        m_minTrkDistY.push_back(minDeltaY);		

        if ((*hitTruthPdgId)[i].size()>0) {
          if ((*hitTruthPdgId)[i][j].size()>0) {
            int nPUinCluster = 0;
            for (int k=0; k<(int)(*hitTruthPdgId)[i][j].size(); k++) {
              if ((*hitTruthPdgId)[i][j][k]>2.0e5) { nPUinCluster++; }
              m_clus_contributingPtcsBarcode.push_back((*hitTruthPdgId)[i][j][k]);
            }
            m_clus_nContributingPtcs.push_back((*hitTruthPdgId)[i][j].size());
            m_clus_nContributingPU.push_back(nPUinCluster);
          }
        }

      }  // end of cluster loop

      m_trktree->Fill();
    }

  }
}

void makeTree::finalize() {
  m_rootFile->Write();
  m_rootFile->Close();
}



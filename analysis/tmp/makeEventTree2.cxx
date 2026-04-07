#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "makeEventTree2.h"
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

void makeEventTree2::initialize(TString rout) {
  m_rootFile = new TFile(rout,"recreate");

  m_eventTree2 = new TTree("eventtree","eventtree");
  m_eventTree2->SetDirectory(m_rootFile);
  m_eventTree2->Branch("runNumber",&m_runNumber);
  m_eventTree2->Branch("eventNumber",&m_eventNumber);
  m_eventTree2->Branch("pileupWeight",&m_pileupWeight);
  m_eventTree2->Branch("mcWeight",&m_mcWeight);
  m_eventTree2->Branch("lumiBlock",&m_lumiBlock);
  m_eventTree2->Branch("averagePU",&m_averagePU);
  m_eventTree2->Branch("eventPU",&m_eventPU);
  m_eventTree2->Branch("mcFlag",&m_mcFlag);
  m_eventTree2->Branch("numVtx",&m_numVtx);
  m_eventTree2->Branch("pVtxX",&m_pVtxX);
  m_eventTree2->Branch("pVtxY",&m_pVtxY);
  m_eventTree2->Branch("pVtxZ",&m_pVtxZ);
  m_eventTree2->Branch("pVtxXErr",&m_pVtxXErr);
  m_eventTree2->Branch("pVtxYErr",&m_pVtxYErr);
  m_eventTree2->Branch("pVtxZErr",&m_pVtxZErr);
  m_eventTree2->Branch("truthPVtxX",&m_truthPVtxX);
  m_eventTree2->Branch("truthPVtxY",&m_truthPVtxY);
  m_eventTree2->Branch("truthPVtxZ",&m_truthPVtxZ);
  m_eventTree2->Branch("trackPt",&m_trackPt);
  m_eventTree2->Branch("trackPhi",&m_trackPhi);
  m_eventTree2->Branch("trackEta",&m_trackEta);
  m_eventTree2->Branch("trackTheta",&m_trackTheta);
  m_eventTree2->Branch("trackCharge",&m_trackCharge);
  m_eventTree2->Branch("trackqOverP",&m_trackqOverP);
  m_eventTree2->Branch("trackD0",&m_trackD0);
  m_eventTree2->Branch("trackZ0",&m_trackZ0);
  m_eventTree2->Branch("trackD0Err",&m_trackD0Err);
  m_eventTree2->Branch("trackZ0Err",&m_trackZ0Err);
  m_eventTree2->Branch("trackqOverPErr",&m_trackqOverPErr);
  m_eventTree2->Branch("trackDeltaZSinTheta",&m_trackDeltaZSinTheta);
  m_eventTree2->Branch("trackMindR",&m_trackMindR);
  m_eventTree2->Branch("trackClass",&m_trackClass);
  m_eventTree2->Branch("trackPassCut",&m_trackPassCut);
  m_eventTree2->Branch("trackOutliers",&m_trackOutliers);
  m_eventTree2->Branch("trackChiSqPerDof",&m_trackChiSqPerDof);
  m_eventTree2->Branch("trackStdDevChi2OS",&m_trackStdDevChi2OS);
  m_eventTree2->Branch("trackNPixelHits",&m_trackNPixelHits);
  m_eventTree2->Branch("trackNPixelHoles",&m_trackNPixelHoles);
  m_eventTree2->Branch("trackNPixelLayers",&m_trackNPixelLayers);
  m_eventTree2->Branch("trackNPixelOutliers",&m_trackNPixelOutliers);
  m_eventTree2->Branch("nPixelGanged",&m_nPixelGanged);
  m_eventTree2->Branch("nPixelShared",&m_nPixelShared);
  m_eventTree2->Branch("nPixelSplit",&m_nPixelSplit);
  m_eventTree2->Branch("nPixelDeadSensors",&m_nPixelDeadSensors);
  m_eventTree2->Branch("trackNIBLUsedHits",&m_trackNIBLUsedHits);
  m_eventTree2->Branch("trackNBLUsedHits",&m_trackNBLUsedHits);
  m_eventTree2->Branch("trackNSCTHits",&m_trackNSCTHits);
  m_eventTree2->Branch("trackNSCTHoles",&m_trackNSCTHoles);
  m_eventTree2->Branch("nSCTOutliers",&m_nSCTOutliers);
  m_eventTree2->Branch("nSCTDeadSensors",&m_nSCTDeadSensors);
  m_eventTree2->Branch("trackNTRTHits",&m_trackNTRTHits);
  m_eventTree2->Branch("nTRTOutliers",&m_nTRTOutliers);
  m_eventTree2->Branch("nTRTHoles",&m_nTRTHoles);
  m_eventTree2->Branch("nTRTHTHits",&m_nTRTHTHits);
  m_eventTree2->Branch("trackNSplitIBLHits",&m_trackNSplitIBLHits);
  m_eventTree2->Branch("trackNSplitBLHits",&m_trackNSplitBLHits);
  m_eventTree2->Branch("trackNSharedIBLHits",&m_trackNSharedIBLHits);
  m_eventTree2->Branch("trackNSharedBLHits",&m_trackNSharedBLHits);
  m_eventTree2->Branch("trackNIBLHits",&m_trackNIBLHits);
  m_eventTree2->Branch("trackNBLHits",&m_trackNBLHits);
  m_eventTree2->Branch("trackNL1Hits",&m_trackNL1Hits);
  m_eventTree2->Branch("trackNL2Hits",&m_trackNL2Hits);
  m_eventTree2->Branch("trackNECHits",&m_trackNECHits);
  m_eventTree2->Branch("trackNIBLExpectedHits",&m_trackNIBLExpectedHits);
  m_eventTree2->Branch("trackNBLExpectedHits",&m_trackNBLExpectedHits);
  m_eventTree2->Branch("trackNSCTSharedHits",&m_trackNSCTSharedHits);
  m_eventTree2->Branch("trackTruncateddEdx",&m_trackTruncateddEdx);
  m_eventTree2->Branch("trackNTruncateddEdx",&m_trackNTruncateddEdx);
  m_eventTree2->Branch("extrapolatedIBLX",&m_extrapolatedIBLX);
  m_eventTree2->Branch("extrapolatedIBLY",&m_extrapolatedIBLY);
  m_eventTree2->Branch("extrapolatedIBLZ",&m_extrapolatedIBLZ);
  m_eventTree2->Branch("extrapolatedBLX",&m_extrapolatedBLX);
  m_eventTree2->Branch("extrapolatedBLY",&m_extrapolatedBLY);
  m_eventTree2->Branch("extrapolatedBLZ",&m_extrapolatedBLZ);
  m_eventTree2->Branch("extrapolatedL1X",&m_extrapolatedL1X);
  m_eventTree2->Branch("extrapolatedL1Y",&m_extrapolatedL1Y);
  m_eventTree2->Branch("extrapolatedL1Z",&m_extrapolatedL1Z);
  m_eventTree2->Branch("extrapolatedL2X",&m_extrapolatedL2X);
  m_eventTree2->Branch("extrapolatedL2Y",&m_extrapolatedL2Y);
  m_eventTree2->Branch("extrapolatedL2Z",&m_extrapolatedL2Z);
  m_eventTree2->Branch("truePt",&m_truePt);
  m_eventTree2->Branch("truePhi",&m_truePhi);
  m_eventTree2->Branch("trueEta",&m_trueEta);
  m_eventTree2->Branch("trueE",&m_trueE);
  m_eventTree2->Branch("trued0",&m_trued0);
  m_eventTree2->Branch("truez0",&m_truez0);
  m_eventTree2->Branch("truephi0",&m_truephi0);
  m_eventTree2->Branch("truetheta",&m_truetheta);
  m_eventTree2->Branch("trueqoverp",&m_trueqoverp);
  m_eventTree2->Branch("truepdgid",&m_truepdgid);
  m_eventTree2->Branch("truebarcode",&m_truebarcode);
  m_eventTree2->Branch("truthmatchprob",&m_truthmatchprob);
  m_eventTree2->Branch("trueMindR",&m_trueMindR);
  m_eventTree2->Branch("genVtxR",&m_genVtxR);
  m_eventTree2->Branch("genVtxZ",&m_genVtxZ);
  m_eventTree2->Branch("parentFlavour",&m_parentFlavour);
  m_eventTree2->Branch("hitIsEndCap",&m_hitIsEndCap);
  m_eventTree2->Branch("hitIsHole",&m_hitIsHole);
  m_eventTree2->Branch("hitIsOutlier",&m_hitIsOutlier);
  m_eventTree2->Branch("hitLayer",&m_hitLayer);
  m_eventTree2->Branch("hitEtaModule",&m_hitEtaModule);
  m_eventTree2->Branch("hitPhiModule",&m_hitPhiModule);
  m_eventTree2->Branch("hitCharge",&m_hitCharge);
  m_eventTree2->Branch("hitToT",&m_hitToT);
  m_eventTree2->Branch("hitLVL1A",&m_hitLVL1A);
  m_eventTree2->Branch("hitNPixel",&m_hitNPixel);
  m_eventTree2->Branch("hitNPixelX",&m_hitNPixelX);
  m_eventTree2->Branch("hitNPixelY",&m_hitNPixelY);
  m_eventTree2->Branch("hitBSerr",&m_hitBSerr);
  m_eventTree2->Branch("hitDCSstate",&m_hitDCSstate);
  m_eventTree2->Branch("hitVBias",&m_hitVBias);
  m_eventTree2->Branch("hitTemp",&m_hitTemp);
  m_eventTree2->Branch("hitLorentzShift",&m_hitLorentzShift);
  m_eventTree2->Branch("hitIsSplit",&m_hitIsSplit);
  m_eventTree2->Branch("hitIsEdge",&m_hitIsEdge);
  m_eventTree2->Branch("hitIsOverflow",&m_hitIsOverflow);
  m_eventTree2->Branch("hitIsolation10x2",&m_hitIsolation10x2);
  m_eventTree2->Branch("hitIsolation20x4",&m_hitIsolation20x4);
  m_eventTree2->Branch("hitGlobalX",&m_hitGlobalX);
  m_eventTree2->Branch("hitGlobalY",&m_hitGlobalY);
  m_eventTree2->Branch("hitGlobalZ",&m_hitGlobalZ);
  m_eventTree2->Branch("hitLocalX",&m_hitLocalX);
  m_eventTree2->Branch("hitLocalY",&m_hitLocalY);
  m_eventTree2->Branch("trkLocalX",&m_trkLocalX);
  m_eventTree2->Branch("trkLocalY",&m_trkLocalY);
  m_eventTree2->Branch("unbiasedResidualX",&m_unbiasedResidualX);
  m_eventTree2->Branch("unbiasedResidualY",&m_unbiasedResidualY);
  m_eventTree2->Branch("unbiasedPullX",&m_unbiasedPullX);
  m_eventTree2->Branch("unbiasedPullY",&m_unbiasedPullY);
  m_eventTree2->Branch("trkPhiOnSurface",&m_trkPhiOnSurface);
  m_eventTree2->Branch("trkThetaOnSurface",&m_trkThetaOnSurface);
  m_eventTree2->Branch("numTotalClustersPerModule",&m_numTotalClustersPerModule);
  m_eventTree2->Branch("numTotalPixelsPerModule",&m_numTotalPixelsPerModule);
  m_eventTree2->Branch("RdoToT",&m_RdoToT);
  m_eventTree2->Branch("RdoCharge",&m_RdoCharge);
  m_eventTree2->Branch("RdoPhi",&m_RdoPhi);
  m_eventTree2->Branch("RdoEta",&m_RdoEta);
  m_eventTree2->Branch("hitTruthPdgId",&m_hitTruthPdgId);
  m_eventTree2->Branch("siHitPdgId",&m_siHitPdgId);
  m_eventTree2->Branch("siHitBarcode",&m_siHitBarcode);
  m_eventTree2->Branch("siHitStartPosX",&m_siHitStartPosX);
  m_eventTree2->Branch("siHitStartPosY",&m_siHitStartPosY);
  m_eventTree2->Branch("siHitEndPosX",&m_siHitEndPosX);
  m_eventTree2->Branch("siHitEndPosY",&m_siHitEndPosY);
  m_eventTree2->Branch("siHitEnergyDeposit",&m_siHitEnergyDeposit);
  m_eventTree2->Branch("truthEta",&m_truthEta);
  m_eventTree2->Branch("truthPhi",&m_truthPhi);
  m_eventTree2->Branch("truthPt",&m_truthPt);
  m_eventTree2->Branch("truthE",&m_truthE);
  m_eventTree2->Branch("truthCharge",&m_truthCharge);
  m_eventTree2->Branch("truthPdgId",&m_truthPdgId);
  m_eventTree2->Branch("truthBarcode",&m_truthBarcode);
  m_eventTree2->Branch("truthStatus",&m_truthStatus);
  m_eventTree2->Branch("muonEta",&m_muonEta);
  m_eventTree2->Branch("muonPhi",&m_muonPhi);
  m_eventTree2->Branch("muonPt",&m_muonPt);
  m_eventTree2->Branch("muonE",&m_muonE);
  m_eventTree2->Branch("jetEta",&m_jetEta);
  m_eventTree2->Branch("jetPhi",&m_jetPhi);
  m_eventTree2->Branch("jetPt",&m_jetPt);
  m_eventTree2->Branch("jetE",&m_jetE);
  m_eventTree2->Branch("jetPassJVT",&m_jetPassJVT);
  m_eventTree2->Branch("tauEta",&m_tauEta);
  m_eventTree2->Branch("tauPhi",&m_tauPhi);
  m_eventTree2->Branch("tauPt",&m_tauPt);
  m_eventTree2->Branch("tauM",&m_tauM);
  m_eventTree2->Branch("tauCharge",&m_tauCharge);
  m_eventTree2->Branch("tauNumTracksCharged",&m_tauNumTracksCharged);
  m_eventTree2->Branch("tauNumTracksIsolation",&m_tauNumTracksIsolation);
  m_eventTree2->Branch("tauNumAllTracks",&m_tauNumAllTracks);
  m_eventTree2->Branch("tauTrackFilterProngs",&m_tauTrackFilterProngs);
  m_eventTree2->Branch("tauTrackFilterQuality",&m_tauTrackFilterQuality);
  m_eventTree2->Branch("tauNumClusters",&m_tauNumClusters);
  m_eventTree2->Branch("tauNumPi0s",&m_tauNumPi0s);
  m_eventTree2->Branch("tauNumCells",&m_tauNumCells);
  m_eventTree2->Branch("tauIsTauFlags",&m_tauIsTauFlags);
  m_eventTree2->Branch("tauBDTJetScore",&m_tauBDTJetScore);
  m_eventTree2->Branch("tauBDTEleScore",&m_tauBDTEleScore);
  m_eventTree2->Branch("tauEleMatchLikelihoodScore",&m_tauEleMatchLikelihoodScore);
  m_eventTree2->Branch("tauBDTJetScoreSigTrans",&m_tauBDTJetScoreSigTrans);
  m_eventTree2->Branch("tauRNNJetScore",&m_tauRNNJetScore);
  m_eventTree2->Branch("tauRNNJetScoreSigTrans",&m_tauRNNJetScoreSigTrans);
  m_eventTree2->Branch("tauIpZ0SinThetaSigLeadTrk",&m_tauIpZ0SinThetaSigLeadTrk);
  m_eventTree2->Branch("tauEtOverPtLeadTrk",&m_tauEtOverPtLeadTrk);
  m_eventTree2->Branch("tauLeadTrkPt",&m_tauLeadTrkPt);
  m_eventTree2->Branch("tauIpSigLeadTrk",&m_tauIpSigLeadTrk);
  m_eventTree2->Branch("tauMassTrkSys",&m_tauMassTrkSys);
  m_eventTree2->Branch("tauTrkWidth2",&m_tauTrkWidth2);
  m_eventTree2->Branch("tauTrFlightPathSig",&m_tauTrFlightPathSig);
  m_eventTree2->Branch("tauEMRadius",&m_tauEMRadius);
  m_eventTree2->Branch("tauHadRadius",&m_tauHadRadius);
  m_eventTree2->Branch("tauEtEMAtEMScale",&m_tauEtEMAtEMScale);
  m_eventTree2->Branch("tauEtHadAtEMScale",&m_tauEtHadAtEMScale);
  m_eventTree2->Branch("tauIsolFrac",&m_tauIsolFrac);
  m_eventTree2->Branch("tauCentFrac",&m_tauCentFrac);
  m_eventTree2->Branch("tauStripWidth2",&m_tauStripWidth2);
  m_eventTree2->Branch("tauNumStrip",&m_tauNumStrip);
  m_eventTree2->Branch("tauTrkAvgDist",&m_tauTrkAvgDist);
  m_eventTree2->Branch("tauDRmax",&m_tauDRmax);
  m_eventTree2->Branch("metCaloPx",&m_metCaloPx);
  m_eventTree2->Branch("metCaloPy",&m_metCaloPy);
  m_eventTree2->Branch("metCaloET",&m_metCaloET);
  m_eventTree2->Branch("metCaloSumET",&m_metCaloSumET);
  m_eventTree2->Branch("metTrackPx",&m_metTrackPx);
  m_eventTree2->Branch("metTrackPy",&m_metTrackPy);
  m_eventTree2->Branch("metTrackET",&m_metTrackET);
  m_eventTree2->Branch("metTrackSumET",&m_metTrackSumET);
  m_eventTree2->Branch("metLocHadPx",&m_metLocHadPx);
  m_eventTree2->Branch("metLocHadPy",&m_metLocHadPy);
  m_eventTree2->Branch("metLocHadET",&m_metLocHadET);
  m_eventTree2->Branch("metLocHadSumET",&m_metLocHadSumET);
  m_eventTree2->Branch("metCoreLCTopoPx",&m_metCoreLCTopoPx);
  m_eventTree2->Branch("metCoreLCTopoPy",&m_metCoreLCTopoPy);
  m_eventTree2->Branch("metCoreLCTopoET",&m_metCoreLCTopoET);
  m_eventTree2->Branch("metCoreLCTopoSumET",&m_metCoreLCTopoSumET);
  m_eventTree2->Branch("metRefLCTopoPx",&m_metRefLCTopoPx);
  m_eventTree2->Branch("metRefLCTopoPy",&m_metRefLCTopoPy);
  m_eventTree2->Branch("metRefLCTopoET",&m_metRefLCTopoET);
  m_eventTree2->Branch("metRefLCTopoSumET",&m_metRefLCTopoSumET);
  m_eventTree2->Branch("metTruthPx",&m_metTruthPx);
  m_eventTree2->Branch("metTruthPy",&m_metTruthPy);
  m_eventTree2->Branch("metTruthET",&m_metTruthET);
  m_eventTree2->Branch("metTruthSumET",&m_metTruthSumET);

}

void makeEventTree2::execute() {

  Long64_t nentries = fChain->GetEntriesFast();

  cout << "Total Number of Evetns = " << nentries << endl;

  Int_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (jentry%10000==0) { cout << jentry << " events processed..." << endl; }

    TLorentzVector muon1,muon2,Zboson;
    muon1.SetPtEtaPhiE((*muonPt)[0],(*muonEta)[0],(*muonPhi)[0],(*muonE)[0]);
    muon2.SetPtEtaPhiE((*muonPt)[1],(*muonEta)[1],(*muonPhi)[1],(*muonE)[1]);
    Zboson = muon1 + muon2;

    m_runNumber    = runNumber;
    m_eventNumber  = eventNumber;
    m_pileupWeight = pileupWeight;
    m_mcWeight     = mcWeight;
    m_lumiBlock    = lumiBlock;
    m_averagePU    = averagePU;
    m_eventPU      = eventPU;
    m_mcFlag       = mcFlag;
    m_numVtx       = numVtx;
    m_pVtxX        = pVtxX;
    m_pVtxY        = pVtxY;
    m_pVtxZ        = pVtxZ;
    m_pVtxXErr     = pVtxXErr;
    m_pVtxYErr     = pVtxYErr;
    m_pVtxZErr     = pVtxZErr;
    m_truthPVtxX   = truthPVtxX;
    m_truthPVtxY   = truthPVtxY;
    m_truthPVtxZ   = truthPVtxZ;

    m_trackPt.clear();
    m_trackPhi.clear();
    m_trackEta.clear();
    m_trackTheta.clear();
    m_trackCharge.clear();
    m_trackqOverP.clear();
    m_trackD0.clear();
    m_trackZ0.clear();
    m_trackD0Err.clear();
    m_trackZ0Err.clear();
    m_trackqOverPErr.clear();
    m_trackDeltaZSinTheta.clear();
    m_trackMindR.clear();
    m_trackClass.clear();
    m_trackPassCut.clear();
    m_trackOutliers.clear();
    m_trackChiSqPerDof.clear();
    m_trackStdDevChi2OS.clear();
    m_trackNPixelHits.clear();
    m_trackNPixelHoles.clear();
    m_trackNPixelLayers.clear();
    m_trackNPixelOutliers.clear();
    m_nPixelGanged.clear();
    m_nPixelShared.clear();
    m_nPixelSplit.clear();
    m_nPixelDeadSensors.clear();
    m_trackNIBLUsedHits.clear();
    m_trackNBLUsedHits.clear();
    m_trackNSCTHits.clear();
    m_trackNSCTHoles.clear();
    m_nSCTOutliers.clear();
    m_nSCTDeadSensors.clear();
    m_trackNTRTHits.clear();
    m_nTRTOutliers.clear();
    m_nTRTHoles.clear();
    m_nTRTHTHits.clear();
    m_trackNSplitIBLHits.clear();
    m_trackNSplitBLHits.clear();
    m_trackNSharedIBLHits.clear();
    m_trackNSharedBLHits.clear();
    m_trackNIBLHits.clear();
    m_trackNBLHits.clear();
    m_trackNL1Hits.clear();
    m_trackNL2Hits.clear();
    m_trackNECHits.clear();
    m_trackNIBLExpectedHits.clear();
    m_trackNBLExpectedHits.clear();
    m_trackNSCTSharedHits.clear();
    m_trackTruncateddEdx.clear();
    m_trackNTruncateddEdx.clear();
    m_extrapolatedIBLX.clear();
    m_extrapolatedIBLY.clear();
    m_extrapolatedIBLZ.clear();
    m_extrapolatedBLX.clear();
    m_extrapolatedBLY.clear();
    m_extrapolatedBLZ.clear();
    m_extrapolatedL1X.clear();
    m_extrapolatedL1Y.clear();
    m_extrapolatedL1Z.clear();
    m_extrapolatedL2X.clear();
    m_extrapolatedL2Y.clear();
    m_extrapolatedL2Z.clear();
    m_truePt.clear();
    m_truePhi.clear();
    m_trueEta.clear();
    m_trueE.clear();
    m_trued0.clear();
    m_truez0.clear();
    m_truephi0.clear();
    m_truetheta.clear();
    m_trueqoverp.clear();
    m_truepdgid.clear();
    m_truebarcode.clear();
    m_truthmatchprob.clear();
    m_trueMindR.clear();
    m_genVtxR.clear();
    m_genVtxZ.clear();
    m_parentFlavour.clear();
    m_hitIsEndCap.clear();
    m_hitIsHole.clear();
    m_hitIsOutlier.clear();
    m_hitLayer.clear();
    m_hitEtaModule.clear();
    m_hitPhiModule.clear();
    m_hitCharge.clear();
    m_hitToT.clear();
    m_hitLVL1A.clear();
    m_hitNPixel.clear();
    m_hitNPixelX.clear();
    m_hitNPixelY.clear();
    m_hitBSerr.clear();
    m_hitDCSstate.clear();
    m_hitVBias.clear();
    m_hitTemp.clear();
    m_hitLorentzShift.clear();
    m_hitIsSplit.clear();
    m_hitIsEdge.clear();
    m_hitIsOverflow.clear();
    m_hitIsolation10x2.clear();
    m_hitIsolation20x4.clear();
    m_hitGlobalX.clear();
    m_hitGlobalY.clear();
    m_hitGlobalZ.clear();
    m_hitLocalX.clear();
    m_hitLocalY.clear();
    m_trkLocalX.clear();
    m_trkLocalY.clear();
    m_unbiasedResidualX.clear();
    m_unbiasedResidualY.clear();
    m_unbiasedPullX.clear();
    m_unbiasedPullY.clear();
    m_trkPhiOnSurface.clear();
    m_trkThetaOnSurface.clear();
    m_numTotalClustersPerModule.clear();
    m_numTotalPixelsPerModule.clear();
    m_RdoToT.clear();
    m_RdoCharge.clear();
    m_RdoPhi.clear();
    m_RdoEta.clear();
    m_hitTruthPdgId.clear();
    m_siHitPdgId.clear();
    m_siHitBarcode.clear();
    m_siHitStartPosX.clear();
    m_siHitStartPosY.clear();
    m_siHitEndPosX.clear();
    m_siHitEndPosY.clear();
    m_siHitEnergyDeposit.clear();
    for (int i=0; i<(int)(*trackPt).size(); i++) {

      TLorentzVector trk1;
      trk1.SetPtEtaPhiM((*trackPt)[i],(*trackEta)[i],(*trackPhi)[i],0.139);
      if (trk1.DeltaR(muon1)>0.8 && trk1.DeltaR(muon2)>0.8) { continue; }


/*
//      if ((*trackPt)[i]<2.0) { continue; }
      // Loose track selection
      bool trackSel = false;
      if ((TMath::Abs((*trackEta)[i])<=1.65 && (*trackNSCTHits)[i]>=6) || (TMath::Abs((*trackEta)[i])>1.65 && (*trackNSCTHits)[i]>=8)) { trackSel=true; }
      if (!trackSel) { continue; }

      // Track quality
      if ((*trackNPixelHits)[i]<4)     { continue; }
      if ((*trackNPixelHoles)[i]>0)    { continue; }
      if ((*nPixelDeadSensors)[i]>0)   { continue; }
      if ((*trackNPixelOutliers)[i]>0) { continue; }
      if ((*nPixelShared)[i]>0)        { continue; }
      if ((*nPixelSplit)[i]>0)         { continue; } 

      // Pixel hit multiplicity
      bool pixelHit = false;
      if      ((*trackNIBLHits)[i]>0 && (*trackNBLHits)[i]>0 && (*trackNL1Hits)[i]>0 && (*trackNL2Hits)[i]>0) { pixelHit=true; }
      else if ((*trackNIBLHits)[i]>0 && (*trackNECHits)[i]>1) { pixelHit=true; }
      if (!pixelHit) { continue; }

      // Cluster should not contain edge pixels
      int isEdge = 0;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) { isEdge+=(*hitIsEdge)[i][j]; }
      if (isEdge>0) { continue; }

      // Cluster should not contain overflow pixels
      int isOverflow = 0;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) { isOverflow+=(*hitIsOverflow)[i][j]; }
      if (isOverflow>0) { continue; }

      // Cluster should not contain split state
      int isSplit = 0;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) { isSplit+=(*hitIsSplit)[i][j]; }
      if (isSplit>0) { continue; }

      // Strong isolation
      int iso20x4 = 0;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) { iso20x4+=(*hitIsolation20x4)[i][j]; }
      if (iso20x4>0) { continue; }

      // Cut on angle alpha
      bool isAlphaCut = false;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {
        float alpha = TMath::ATan(TMath::Sqrt(TMath::Power(TMath::Tan((*trkThetaOnSurface)[i][j]),2)+TMath::Power(TMath::Tan((*trkPhiOnSurface)[i][j]),2)));
        if (TMath::Cos(alpha)<0.16) { isAlphaCut=true; break; }
      }
      if (isAlphaCut) { continue; }

      // Exclude 3D sensor since no Lorentz angle, different thickness (230um) than IBL
      bool exclude3D = false;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {
        if ((*hitIsEndCap)[i][j]==0 && (*hitLayer)[i][j]==0) {  // IBL
          if ((*hitEtaModule)[i][j]>5)  { exclude3D=true; break; }
          if ((*hitEtaModule)[i][j]<-6) { exclude3D=true; break; }
        }
      }
      if (exclude3D) { continue; }

      // Exclude small hit ToT<10 for pixel (not IBL)
      bool excludeSmallToT = false;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {
        if ((*hitIsEndCap)[i][j]==0 && (*hitLayer)[i][j]>0) {  // Barrel NOT IBL
          if ((*hitToT)[i][j]<10) { excludeSmallToT=true; break; }
        }
      }
      if (excludeSmallToT) { continue; }
*/



      m_trackPt.push_back((*trackPt)[i]);
      m_trackPhi.push_back((*trackPhi)[i]);
      m_trackEta.push_back((*trackEta)[i]);
      m_trackTheta.push_back((*trackTheta)[i]);
      m_trackCharge.push_back((*trackCharge)[i]);
      m_trackqOverP.push_back((*trackqOverP)[i]);
      m_trackD0.push_back((*trackD0)[i]);
      m_trackZ0.push_back((*trackZ0)[i]);
      m_trackD0Err.push_back((*trackD0Err)[i]);
      m_trackZ0Err.push_back((*trackZ0Err)[i]);
      m_trackqOverPErr.push_back((*trackqOverPErr)[i]);
      m_trackDeltaZSinTheta.push_back((*trackDeltaZSinTheta)[i]);
      m_trackMindR.push_back((*trackMindR)[i]);
      m_trackClass.push_back((*trackClass)[i]);
      m_trackPassCut.push_back((*trackPassCut)[i]);
      m_trackOutliers.push_back((*trackOutliers)[i]);
      m_trackChiSqPerDof.push_back((*trackChiSqPerDof)[i]);
      m_trackStdDevChi2OS.push_back((*trackStdDevChi2OS)[i]);
      m_trackNPixelHits.push_back((*trackNPixelHits)[i]);
      m_trackNPixelHoles.push_back((*trackNPixelHoles)[i]);
      m_trackNPixelLayers.push_back((*trackNPixelLayers)[i]);
      m_trackNPixelOutliers.push_back((*trackNPixelOutliers)[i]);
      m_nPixelGanged.push_back((*nPixelGanged)[i]);
      m_nPixelShared.push_back((*nPixelShared)[i]);
      m_nPixelSplit.push_back((*nPixelSplit)[i]);
      m_nPixelDeadSensors.push_back((*nPixelDeadSensors)[i]);
      m_trackNIBLUsedHits.push_back((*trackNIBLUsedHits)[i]);
      m_trackNBLUsedHits.push_back((*trackNBLUsedHits)[i]);
      m_trackNSCTHits.push_back((*trackNSCTHits)[i]);
      m_trackNSCTHoles.push_back((*trackNSCTHoles)[i]);
      m_nSCTOutliers.push_back((*nSCTOutliers)[i]);
      m_nSCTDeadSensors.push_back((*nSCTDeadSensors)[i]);
      m_trackNTRTHits.push_back((*trackNTRTHits)[i]);
      m_nTRTOutliers.push_back((*nTRTOutliers)[i]);
      m_nTRTHoles.push_back((*nTRTHoles)[i]);
      m_nTRTHTHits.push_back((*nTRTHTHits)[i]);
      m_trackNSplitIBLHits.push_back((*trackNSplitIBLHits)[i]);
      m_trackNSplitBLHits.push_back((*trackNSplitBLHits)[i]);
      m_trackNSharedIBLHits.push_back((*trackNSharedIBLHits)[i]);
      m_trackNSharedBLHits.push_back((*trackNSharedBLHits)[i]);
      m_trackNIBLHits.push_back((*trackNIBLHits)[i]);
      m_trackNBLHits.push_back((*trackNBLHits)[i]);
      m_trackNL1Hits.push_back((*trackNL1Hits)[i]);
      m_trackNL2Hits.push_back((*trackNL2Hits)[i]);
      m_trackNECHits.push_back((*trackNECHits)[i]);
      m_trackNIBLExpectedHits.push_back((*trackNIBLExpectedHits)[i]);
      m_trackNBLExpectedHits.push_back((*trackNBLExpectedHits)[i]);
      m_trackNSCTSharedHits.push_back((*trackNSCTSharedHits)[i]);
      m_trackTruncateddEdx.push_back((*trackTruncateddEdx)[i]);
      m_trackNTruncateddEdx.push_back((*trackNTruncateddEdx)[i]);
      m_extrapolatedIBLX.push_back((*extrapolatedIBLX)[i]);
      m_extrapolatedIBLY.push_back((*extrapolatedIBLY)[i]);
      m_extrapolatedIBLZ.push_back((*extrapolatedIBLZ)[i]);
      m_extrapolatedBLX.push_back((*extrapolatedBLX)[i]);
      m_extrapolatedBLY.push_back((*extrapolatedBLY)[i]);
      m_extrapolatedBLZ.push_back((*extrapolatedBLZ)[i]);
      m_extrapolatedL1X.push_back((*extrapolatedL1X)[i]);
      m_extrapolatedL1Y.push_back((*extrapolatedL1Y)[i]);
      m_extrapolatedL1Z.push_back((*extrapolatedL1Z)[i]);
      m_extrapolatedL2X.push_back((*extrapolatedL2X)[i]);
      m_extrapolatedL2Y.push_back((*extrapolatedL2Y)[i]);
      m_extrapolatedL2Z.push_back((*extrapolatedL2Z)[i]);
      m_truePt.push_back((*truePt)[i]);
      m_truePhi.push_back((*truePhi)[i]);
      m_trueEta.push_back((*trueEta)[i]);
      m_trueE.push_back((*trueE)[i]);
      m_trued0.push_back((*trued0)[i]);
      m_truez0.push_back((*truez0)[i]);
      m_truephi0.push_back((*truephi0)[i]);
      m_truetheta.push_back((*truetheta)[i]);
      m_trueqoverp.push_back((*trueqoverp)[i]);
      m_truepdgid.push_back((*truepdgid)[i]);
      m_truebarcode.push_back((*truebarcode)[i]);
      m_truthmatchprob.push_back((*truthmatchprob)[i]);
      m_trueMindR.push_back((*trueMindR)[i]);
      m_genVtxR.push_back((*genVtxR)[i]);
      m_genVtxZ.push_back((*genVtxZ)[i]);
      m_parentFlavour.push_back((*parentFlavour)[i]);

      vector<int> hitIsEndCapTmp;
      vector<int> hitIsHoleTmp;
      vector<int> hitIsOutlierTmp;
      vector<int> hitLayerTmp;
      vector<int> hitEtaModuleTmp;
      vector<int> hitPhiModuleTmp;
      vector<float> hitChargeTmp;
      vector<int> hitToTTmp;
      vector<int> hitLVL1ATmp;
      vector<int> hitNPixelTmp;
      vector<int> hitNPixelXTmp;
      vector<int> hitNPixelYTmp;
      vector<int> hitBSerrTmp;
      vector<int> hitDCSstateTmp;
      vector<float> hitVBiasTmp;
      vector<float> hitTempTmp;
      vector<float> hitLorentzShiftTmp;
      vector<int> hitIsSplitTmp;
      vector<bool> hitIsEdgeTmp;
      vector<bool> hitIsOverflowTmp;
      vector<int> hitIsolation10x2Tmp;
      vector<int> hitIsolation20x4Tmp;
      vector<float> hitGlobalXTmp;
      vector<float> hitGlobalYTmp;
      vector<float> hitGlobalZTmp;
      vector<float> hitLocalXTmp;
      vector<float> hitLocalYTmp;
      vector<float> trkLocalXTmp;
      vector<float> trkLocalYTmp;
      vector<float> unbiasedResidualXTmp;
      vector<float> unbiasedResidualYTmp;
      vector<float> unbiasedPullXTmp;
      vector<float> unbiasedPullYTmp;
      vector<float> trkPhiOnSurfaceTmp;
      vector<float> trkThetaOnSurfaceTmp;
      vector<int> numTotalClustersPerModuleTmp;
      vector<int> numTotalPixelsPerModuleTmp;
      for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {
        hitIsEndCapTmp.push_back((*hitIsEndCap)[i][j]);
        hitLayerTmp.push_back((*hitLayer)[i][j]);
        hitEtaModuleTmp.push_back((*hitEtaModule)[i][j]);
        hitPhiModuleTmp.push_back((*hitPhiModule)[i][j]);
        hitChargeTmp.push_back((*hitCharge)[i][j]);
        hitToTTmp.push_back((*hitToT)[i][j]);
        hitLVL1ATmp.push_back((*hitLVL1A)[i][j]);
        hitNPixelTmp.push_back((*hitNPixel)[i][j]);
        hitNPixelXTmp.push_back((*hitNPixelX)[i][j]);
        hitNPixelYTmp.push_back((*hitNPixelY)[i][j]);
        if ((*hitBSerr)[i].size()==(*hitLayer)[i].size()) { hitBSerrTmp.push_back((*hitBSerr)[i][j]); }
        else                                              { hitBSerrTmp.push_back(1); }
        if ((*hitDCSstate)[i].size()==(*hitLayer)[i].size()) { hitDCSstateTmp.push_back((*hitDCSstate)[i][j]); }
        else                                                 { hitDCSstateTmp.push_back(0); }
        hitVBiasTmp.push_back((*hitVBias)[i][j]);
        hitTempTmp.push_back((*hitTemp)[i][j]);
        hitLorentzShiftTmp.push_back((*hitLorentzShift)[i][j]);
        hitIsSplitTmp.push_back((*hitIsSplit)[i][j]);
        hitIsEdgeTmp.push_back((*hitIsEdge)[i][j]);
        hitIsOverflowTmp.push_back((*hitIsOverflow)[i][j]);
        hitIsolation10x2Tmp.push_back((*hitIsolation10x2)[i][j]);
        hitIsolation20x4Tmp.push_back((*hitIsolation20x4)[i][j]);
        hitGlobalXTmp.push_back((*hitGlobalX)[i][j]);
        hitGlobalYTmp.push_back((*hitGlobalY)[i][j]);
        hitGlobalZTmp.push_back((*hitGlobalZ)[i][j]);
        hitLocalXTmp.push_back((*hitLocalX)[i][j]);
        hitLocalYTmp.push_back((*hitLocalY)[i][j]);
        trkLocalXTmp.push_back((*trkLocalX)[i][j]);
        trkLocalYTmp.push_back((*trkLocalY)[i][j]);
        unbiasedResidualXTmp.push_back((*unbiasedResidualX)[i][j]);
        unbiasedResidualYTmp.push_back((*unbiasedResidualY)[i][j]);
        unbiasedPullXTmp.push_back((*unbiasedPullX)[i][j]);
        unbiasedPullYTmp.push_back((*unbiasedPullY)[i][j]);
        trkPhiOnSurfaceTmp.push_back((*trkPhiOnSurface)[i][j]);
        trkThetaOnSurfaceTmp.push_back((*trkThetaOnSurface)[i][j]);
        numTotalClustersPerModuleTmp.push_back((*numTotalClustersPerModule)[i][j]);
        numTotalPixelsPerModuleTmp.push_back((*numTotalPixelsPerModule)[i][j]);
      }

      for (int j=0; j<(int)(*hitIsHole)[i].size(); j++) {
        hitIsHoleTmp.push_back((*hitIsHole)[i][j]);
      }

      for (int j=0; j<(int)(*hitIsOutlier)[i].size(); j++) {
        hitIsOutlierTmp.push_back((*hitIsOutlier)[i][j]);
      }

      m_hitIsEndCap.push_back(hitIsEndCapTmp);
      m_hitIsHole.push_back(hitIsHoleTmp);
      m_hitIsOutlier.push_back(hitIsOutlierTmp);
      m_hitLayer.push_back(hitLayerTmp);
      m_hitEtaModule.push_back(hitEtaModuleTmp);
      m_hitPhiModule.push_back(hitPhiModuleTmp);
      m_hitCharge.push_back(hitChargeTmp);
      m_hitToT.push_back(hitToTTmp);
      m_hitLVL1A.push_back(hitLVL1ATmp);
      m_hitNPixel.push_back(hitNPixelTmp);
      m_hitNPixelX.push_back(hitNPixelXTmp);
      m_hitNPixelY.push_back(hitNPixelYTmp);
      m_hitBSerr.push_back(hitBSerrTmp);
      m_hitDCSstate.push_back(hitDCSstateTmp);
      m_hitVBias.push_back(hitVBiasTmp);
      m_hitTemp.push_back(hitTempTmp);
      m_hitLorentzShift.push_back(hitLorentzShiftTmp);
      m_hitIsSplit.push_back(hitIsSplitTmp);
      m_hitIsEdge.push_back(hitIsEdgeTmp);
      m_hitIsOverflow.push_back(hitIsOverflowTmp);
      m_hitIsolation10x2.push_back(hitIsolation10x2Tmp);
      m_hitIsolation20x4.push_back(hitIsolation20x4Tmp);
      m_hitGlobalX.push_back(hitGlobalXTmp);
      m_hitGlobalY.push_back(hitGlobalYTmp);
      m_hitGlobalZ.push_back(hitGlobalZTmp);
      m_hitLocalX.push_back(hitLocalXTmp);
      m_hitLocalY.push_back(hitLocalYTmp);
      m_trkLocalX.push_back(trkLocalXTmp);
      m_trkLocalY.push_back(trkLocalYTmp);
      m_unbiasedResidualX.push_back(unbiasedResidualXTmp);
      m_unbiasedResidualY.push_back(unbiasedResidualYTmp);
      m_unbiasedPullX.push_back(unbiasedPullXTmp);
      m_unbiasedPullY.push_back(unbiasedPullYTmp);
      m_trkPhiOnSurface.push_back(trkPhiOnSurfaceTmp);
      m_trkThetaOnSurface.push_back(trkThetaOnSurfaceTmp);
      m_numTotalClustersPerModule.push_back(numTotalClustersPerModuleTmp);
      m_numTotalPixelsPerModule.push_back(numTotalPixelsPerModuleTmp);

      vector<vector<int>>   RdoToTTmp2;
      vector<vector<float>> RdoChargeTmp2;
      vector<vector<int>>   RdoPhiTmp2;
      vector<vector<int>>   RdoEtaTmp2;
      vector<vector<int>>   hitTruthPdgIdTmp2;
      if ((*trackPt)[i]>2.0) {
        for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {
          vector<int>   RdoToTTmp;
          vector<float> RdoChargeTmp;
          vector<int>   RdoPhiTmp;
          vector<int>   RdoEtaTmp;
          vector<int>   hitTruthPdgIdTmp;
          for (int k=0; k<(int)(*RdoCharge)[i][j].size(); k++) {
            RdoToTTmp.push_back((*RdoToT)[i][j][k]);
            RdoChargeTmp.push_back((*RdoCharge)[i][j][k]);
            RdoPhiTmp.push_back((*RdoPhi)[i][j][k]);
            RdoEtaTmp.push_back((*RdoEta)[i][j][k]);
            if (mcFlag) {
//              hitTruthPdgIdTmp.push_back((*hitTruthPdgId)[i][j][k]);
            }
          }
          RdoToTTmp2.push_back(RdoToTTmp);
          RdoChargeTmp2.push_back(RdoChargeTmp);
          RdoPhiTmp2.push_back(RdoPhiTmp);
          RdoEtaTmp2.push_back(RdoEtaTmp);
          hitTruthPdgIdTmp2.push_back(hitTruthPdgIdTmp);
        }
      }
      m_RdoToT.push_back(RdoToTTmp2);
      m_RdoCharge.push_back(RdoChargeTmp2);
      m_RdoPhi.push_back(RdoPhiTmp2);
      m_RdoEta.push_back(RdoEtaTmp2);
      m_hitTruthPdgId.push_back(hitTruthPdgIdTmp2);

      vector<vector<int>> siHitPdgIdTmp2;
      vector<vector<int>> siHitBarcodeTmp2;
      vector<vector<float>> siHitStartPosXTmp2;
      vector<vector<float>> siHitStartPosYTmp2;
      vector<vector<float>> siHitEndPosXTmp2;
      vector<vector<float>> siHitEndPosYTmp2;
      vector<vector<float>> siHitEnergyDepositTmp2;
      if ((*siHitPdgId)[i].size()>0 && (*truepdgid)[i]!=0) {
        for (int j=0; j<(int)(*hitLayer)[i].size(); j++) {
          vector<int> siHitPdgIdTmp;
          vector<int> siHitBarcodeTmp;
          vector<float> siHitStartPosXTmp;
          vector<float> siHitStartPosYTmp;
          vector<float> siHitEndPosXTmp;
          vector<float> siHitEndPosYTmp;
          vector<float> siHitEnergyDepositTmp;
          for (int k=0; k<(int)(*siHitPdgId)[i][j].size(); k++) {
            siHitPdgIdTmp.push_back((*siHitPdgId)[i][j][k]);
            siHitBarcodeTmp.push_back((*siHitBarcode)[i][j][k]);
            siHitStartPosXTmp.push_back((*siHitStartPosX)[i][j][k]);
            siHitStartPosYTmp.push_back((*siHitStartPosY)[i][j][k]);
            siHitEndPosXTmp.push_back((*siHitEndPosX)[i][j][k]);
            siHitEndPosYTmp.push_back((*siHitEndPosY)[i][j][k]);
            siHitEnergyDepositTmp.push_back((*siHitEnergyDeposit)[i][j][k]);
          }
          siHitPdgIdTmp2.push_back(siHitPdgIdTmp);
          siHitBarcodeTmp2.push_back(siHitBarcodeTmp);
          siHitStartPosXTmp2.push_back(siHitStartPosXTmp);
          siHitStartPosYTmp2.push_back(siHitStartPosYTmp);
          siHitEndPosXTmp2.push_back(siHitEndPosXTmp);
          siHitEndPosYTmp2.push_back(siHitEndPosYTmp);
          siHitEnergyDepositTmp2.push_back(siHitEnergyDepositTmp);
        }
      }
      m_siHitPdgId.push_back(siHitPdgIdTmp2);
      m_siHitBarcode.push_back(siHitBarcodeTmp2);
      m_siHitStartPosX.push_back(siHitStartPosXTmp2);
      m_siHitStartPosY.push_back(siHitStartPosYTmp2);
      m_siHitEndPosX.push_back(siHitEndPosXTmp2);
      m_siHitEndPosY.push_back(siHitEndPosYTmp2);
      m_siHitEnergyDeposit.push_back(siHitEnergyDepositTmp2);
    }

    m_truthEta.clear();
    m_truthPhi.clear();
    m_truthPt.clear();
    m_truthE.clear();
    m_truthCharge.clear();
    m_truthPdgId.clear();
    m_truthBarcode.clear();
    m_truthStatus.clear();
    for (int i=0; i<(int)(*truthEta).size(); i++) {
      m_truthEta.push_back((*truthEta)[i]);
      m_truthPhi.push_back((*truthPhi)[i]);
      m_truthPt.push_back((*truthPt)[i]);
      m_truthE.push_back((*truthE)[i]);
      m_truthCharge.push_back((*truthCharge)[i]);
      m_truthPdgId.push_back((*truthPdgId)[i]);
      m_truthBarcode.push_back((*truthBarcode)[i]);
      m_truthStatus.push_back((*truthStatus)[i]);
    }

    m_muonEta.clear();
    m_muonPhi.clear();
    m_muonPt.clear();
    m_muonE.clear();
    for (int i=0; i<(int)(*muonPt).size(); i++) {
      m_muonEta.push_back((*muonEta)[i]);
      m_muonPhi.push_back((*muonPhi)[i]);
      m_muonPt.push_back((*muonPt)[i]);
      m_muonE.push_back((*muonE)[i]);
    }

    m_jetEta.clear();
    m_jetPhi.clear();
    m_jetPt.clear();
    m_jetE.clear();
    m_jetPassJVT.clear();
    for (int i=0; i<(int)(*jetPt).size(); i++) {
      m_jetEta.push_back((*jetEta)[i]);
      m_jetPhi.push_back((*jetPhi)[i]);
      m_jetPt.push_back((*jetPt)[i]);
      m_jetE.push_back((*jetE)[i]);
      m_jetPassJVT.push_back((*jetPassJVT)[i]);
    }

    m_tauEta.clear();
    m_tauPhi.clear();
    m_tauPt.clear();
    m_tauM.clear();
    m_tauCharge.clear();
    m_tauNumTracksCharged.clear();
    m_tauNumTracksIsolation.clear();
    m_tauNumAllTracks.clear();
    m_tauTrackFilterProngs.clear();
    m_tauTrackFilterQuality.clear();
    m_tauNumClusters.clear();
    m_tauNumPi0s.clear();
    m_tauNumCells.clear();
    m_tauIsTauFlags.clear();
    m_tauBDTJetScore.clear();
    m_tauBDTEleScore.clear();
    m_tauEleMatchLikelihoodScore.clear();
    m_tauBDTJetScoreSigTrans.clear();
    m_tauRNNJetScore.clear();
    m_tauRNNJetScoreSigTrans.clear();
    m_tauIpZ0SinThetaSigLeadTrk.clear();
    m_tauEtOverPtLeadTrk.clear();
    m_tauLeadTrkPt.clear();
    m_tauIpSigLeadTrk.clear();
    m_tauMassTrkSys.clear();
    m_tauTrkWidth2.clear();
    m_tauTrFlightPathSig.clear();
    m_tauEMRadius.clear();
    m_tauHadRadius.clear();
    m_tauEtEMAtEMScale.clear();
    m_tauEtHadAtEMScale.clear();
    m_tauIsolFrac.clear();
    m_tauCentFrac.clear();
    m_tauStripWidth2.clear();
    m_tauNumStrip.clear();
    m_tauTrkAvgDist.clear();
    m_tauDRmax.clear();
    for (int i=0; i<(int)(*tauPt).size(); i++) {
      m_tauEta.push_back((*tauEta)[i]);
      m_tauPhi.push_back((*tauPhi)[i]);
      m_tauPt.push_back((*tauPt)[i]);
      m_tauM.push_back((*tauM)[i]);
      m_tauCharge.push_back((*tauCharge)[i]);
      m_tauNumTracksCharged.push_back((*tauNumTracksCharged)[i]);
      m_tauNumTracksIsolation.push_back((*tauNumTracksIsolation)[i]);
      m_tauNumAllTracks.push_back((*tauNumAllTracks)[i]);
      m_tauTrackFilterProngs.push_back((*tauTrackFilterProngs)[i]);
      m_tauTrackFilterQuality.push_back((*tauTrackFilterQuality)[i]);
      m_tauNumClusters.push_back((*tauNumClusters)[i]);
      m_tauNumPi0s.push_back((*tauNumPi0s)[i]);
      m_tauNumCells.push_back((*tauNumCells)[i]);
      m_tauIsTauFlags.push_back((*tauIsTauFlags)[i]);
      m_tauBDTJetScore.push_back((*tauBDTJetScore)[i]);
      m_tauBDTEleScore.push_back((*tauBDTEleScore)[i]);
      m_tauEleMatchLikelihoodScore.push_back((*tauEleMatchLikelihoodScore)[i]);
      m_tauBDTJetScoreSigTrans.push_back((*tauBDTJetScoreSigTrans)[i]);
      m_tauRNNJetScore.push_back((*tauRNNJetScore)[i]);
      m_tauRNNJetScoreSigTrans.push_back((*tauRNNJetScoreSigTrans)[i]);
      m_tauIpZ0SinThetaSigLeadTrk.push_back((*tauIpZ0SinThetaSigLeadTrk)[i]);
      m_tauEtOverPtLeadTrk.push_back((*tauEtOverPtLeadTrk)[i]);
      m_tauLeadTrkPt.push_back((*tauLeadTrkPt)[i]);
      m_tauIpSigLeadTrk.push_back((*tauIpSigLeadTrk)[i]);
      m_tauMassTrkSys.push_back((*tauMassTrkSys)[i]);
      m_tauTrkWidth2.push_back((*tauTrkWidth2)[i]);
      m_tauTrFlightPathSig.push_back((*tauTrFlightPathSig)[i]);
      m_tauEMRadius.push_back((*tauEMRadius)[i]);
      m_tauHadRadius.push_back((*tauHadRadius)[i]);
      m_tauEtEMAtEMScale.push_back((*tauEtEMAtEMScale)[i]);
      m_tauEtHadAtEMScale.push_back((*tauEtHadAtEMScale)[i]);
      m_tauIsolFrac.push_back((*tauIsolFrac)[i]);
      m_tauCentFrac.push_back((*tauCentFrac)[i]);
      m_tauStripWidth2.push_back((*tauStripWidth2)[i]);
      m_tauNumStrip.push_back((*tauNumStrip)[i]);
      m_tauTrkAvgDist.push_back((*tauTrkAvgDist)[i]);
      m_tauDRmax.push_back((*tauDRmax)[i]);
    }

    m_metCaloPx          = metCaloPx;
    m_metCaloPy          = metCaloPy;
    m_metCaloET          = metCaloET;
    m_metCaloSumET       = metCaloSumET;
    m_metTrackPx         = metTrackPx;
    m_metTrackPy         = metTrackPy;
    m_metTrackET         = metTrackET;
    m_metTrackSumET      = metTrackSumET;
    m_metLocHadPx        = metLocHadPx;
    m_metLocHadPy        = metLocHadPy;
    m_metLocHadET        = metLocHadET;
    m_metLocHadSumET     = metLocHadSumET;
    m_metCoreLCTopoPx    = metCoreLCTopoPx;
    m_metCoreLCTopoPy    = metCoreLCTopoPy;
    m_metCoreLCTopoET    = metCoreLCTopoET;
    m_metCoreLCTopoSumET = metCoreLCTopoSumET;
    m_metRefLCTopoPx     = metRefLCTopoPx;
    m_metRefLCTopoPy     = metRefLCTopoPy;
    m_metRefLCTopoET     = metRefLCTopoET;
    m_metRefLCTopoSumET  = metRefLCTopoSumET;
    m_metTruthPx         = metTruthPx;
    m_metTruthPy         = metTruthPy;
    m_metTruthET         = metTruthET;
    m_metTruthSumET      = metTruthSumET;

    m_eventTree2->Fill();
  }
}

void makeEventTree2::finalize() {
  m_rootFile->Write();
  m_rootFile->Close();
}


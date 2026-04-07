#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#ifndef eventtree_h
#define eventtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class eventtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   Double_t        pileupWeight;
   Double_t        mcWeight;
   Int_t           lumiBlock;
   Float_t         averagePU;
   Float_t         eventPU;
   Int_t           mcFlag;
   Int_t           numVtx;
   Float_t         pVtxX;
   Float_t         pVtxY;
   Float_t         pVtxZ;
   Float_t         pVtxXErr;
   Float_t         pVtxYErr;
   Float_t         pVtxZErr;
   Float_t         truthPVtxX;
   Float_t         truthPVtxY;
   Float_t         truthPVtxZ;
   vector<float>   *trackPt;
   vector<float>   *trackPhi;
   vector<float>   *trackEta;
   vector<float>   *trackTheta;
   vector<float>   *trackCharge;
   vector<float>   *trackqOverP;
   vector<float>   *trackD0;
   vector<float>   *trackZ0;
   vector<float>   *trackD0Err;
   vector<float>   *trackZ0Err;
   vector<float>   *trackqOverPErr;
   vector<float>   *trackDeltaZSinTheta;
   vector<float>   *trackMindR;
   vector<int>     *trackClass;
   vector<int>     *trackPassCut;
   vector<int>     *trackOutliers;
   vector<float>   *trackChiSqPerDof;
   vector<float>   *trackStdDevChi2OS;
   vector<int>     *trackNPixelHits;
   vector<int>     *trackNPixelHoles;
   vector<int>     *trackNPixelLayers;
   vector<int>     *trackNPixelOutliers;
   vector<int>     *nPixelGanged;
   vector<int>     *nPixelShared;
   vector<int>     *nPixelSplit;
   vector<int>     *nPixelDeadSensors;
   vector<int>     *trackNIBLUsedHits;
   vector<int>     *trackNBLUsedHits;
   vector<int>     *trackNSCTHits;
   vector<int>     *trackNSCTHoles;
   vector<int>     *nSCTOutliers;
   vector<int>     *nSCTDeadSensors;
   vector<int>     *trackNTRTHits;
   vector<int>     *nTRTOutliers;
   vector<int>     *nTRTHoles;
   vector<int>     *nTRTHTHits;
   vector<int>     *trackNSplitIBLHits;
   vector<int>     *trackNSplitBLHits;
   vector<int>     *trackNSharedIBLHits;
   vector<int>     *trackNSharedBLHits;
   vector<int>     *trackNIBLHits;
   vector<int>     *trackNBLHits;
   vector<int>     *trackNL1Hits;
   vector<int>     *trackNL2Hits;
   vector<int>     *trackNECHits;
   vector<int>     *trackNIBLExpectedHits;
   vector<int>     *trackNBLExpectedHits;
   vector<int>     *trackNSCTSharedHits;
   vector<float>   *trackTruncateddEdx;
   vector<int>     *trackNTruncateddEdx;
   vector<float>   *extrapolatedIBLX;
   vector<float>   *extrapolatedIBLY;
   vector<float>   *extrapolatedIBLZ;
   vector<float>   *extrapolatedBLX;
   vector<float>   *extrapolatedBLY;
   vector<float>   *extrapolatedBLZ;
   vector<float>   *extrapolatedL1X;
   vector<float>   *extrapolatedL1Y;
   vector<float>   *extrapolatedL1Z;
   vector<float>   *extrapolatedL2X;
   vector<float>   *extrapolatedL2Y;
   vector<float>   *extrapolatedL2Z;
   vector<float>   *truePt;
   vector<float>   *truePhi;
   vector<float>   *trueEta;
   vector<float>   *trueE;
   vector<float>   *trued0;
   vector<float>   *truez0;
   vector<float>   *truephi0;
   vector<float>   *truetheta;
   vector<float>   *trueqoverp;
   vector<int>     *truepdgid;
   vector<int>     *truebarcode;
   vector<float>   *truthmatchprob;
   vector<float>   *trueMindR;
   vector<float>   *genVtxR;
   vector<float>   *genVtxZ;
   vector<int>     *parentFlavour;
   vector<vector<int> > *hitIsEndCap;
   vector<vector<int> > *hitIsHole;
   vector<vector<int> > *hitIsOutlier;
   vector<vector<int> > *hitLayer;
   vector<vector<int> > *hitEtaModule;
   vector<vector<int> > *hitPhiModule;
   vector<vector<float> > *hitCharge;
   vector<vector<int> > *hitToT;
   vector<vector<int> > *hitLVL1A;
   vector<vector<int> > *hitNPixel;
   vector<vector<int> > *hitNPixelX;
   vector<vector<int> > *hitNPixelY;
   vector<vector<int> > *hitBSerr;
   vector<vector<int> > *hitDCSstate;
   vector<vector<float> > *hitVBias;
   vector<vector<float> > *hitTemp;
   vector<vector<float> > *hitLorentzShift;
   vector<vector<int> > *hitIsSplit;
   vector<vector<bool> > *hitIsEdge;
   vector<vector<bool> > *hitIsOverflow;
   vector<vector<int> > *hitIsolation10x2;
   vector<vector<int> > *hitIsolation20x4;
   vector<vector<float> > *hitGlobalX;
   vector<vector<float> > *hitGlobalY;
   vector<vector<float> > *hitGlobalZ;
   vector<vector<float> > *hitLocalX;
   vector<vector<float> > *hitLocalY;
   vector<vector<float> > *trkLocalX;
   vector<vector<float> > *trkLocalY;
   vector<vector<float> > *unbiasedResidualX;
   vector<vector<float> > *unbiasedResidualY;
   vector<vector<float> > *unbiasedPullX;
   vector<vector<float> > *unbiasedPullY;
   vector<vector<float> > *trkPhiOnSurface;
   vector<vector<float> > *trkThetaOnSurface;
   vector<vector<int> > *numTotalClustersPerModule;
   vector<vector<int> > *numTotalPixelsPerModule;
   vector<vector<vector<int> > > *RdoToT;
   vector<vector<vector<float> > > *RdoCharge;
   vector<vector<vector<int> > > *RdoPhi;
   vector<vector<vector<int> > > *RdoEta;
   vector<vector<vector<int> > > *hitTruthPdgId;
   vector<vector<vector<int> > > *siHitPdgId;
   vector<vector<vector<int> > > *siHitBarcode;
   vector<vector<vector<float> > > *siHitStartPosX;
   vector<vector<vector<float> > > *siHitStartPosY;
   vector<vector<vector<float> > > *siHitEndPosX;
   vector<vector<vector<float> > > *siHitEndPosY;
   vector<vector<vector<float> > > *siHitEnergyDeposit;
   vector<float>   *truthEta;
   vector<float>   *truthPhi;
   vector<float>   *truthPt;
   vector<float>   *truthE;
   vector<float>   *truthCharge;
   vector<int>     *truthPdgId;
   vector<int>     *truthBarcode;
   vector<int>     *truthStatus;
   vector<float>   *muonEta;
   vector<float>   *muonPhi;
   vector<float>   *muonPt;
   vector<float>   *muonE;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetPt;
   vector<float>   *jetE;
   vector<int>     *jetPassJVT;
   vector<float>   *tauEta;
   vector<float>   *tauPhi;
   vector<float>   *tauPt;
   vector<float>   *tauM;
   vector<float>   *tauCharge;
   vector<int>     *tauNumTracksCharged;
   vector<int>     *tauNumTracksIsolation;
   vector<int>     *tauNumAllTracks;
   vector<int>     *tauTrackFilterProngs;
   vector<int>     *tauTrackFilterQuality;
   vector<int>     *tauNumClusters;
   vector<int>     *tauNumPi0s;
   vector<int>     *tauNumCells;
   vector<int>     *tauIsTauFlags;
   vector<float>   *tauBDTJetScore;
   vector<float>   *tauBDTEleScore;
   vector<float>   *tauEleMatchLikelihoodScore;
   vector<float>   *tauBDTJetScoreSigTrans;
   vector<float>   *tauRNNJetScore;
   vector<float>   *tauRNNJetScoreSigTrans;
   vector<float>   *tauIpZ0SinThetaSigLeadTrk;
   vector<float>   *tauEtOverPtLeadTrk;
   vector<float>   *tauLeadTrkPt;
   vector<float>   *tauIpSigLeadTrk;
   vector<float>   *tauMassTrkSys;
   vector<float>   *tauTrkWidth2;
   vector<float>   *tauTrFlightPathSig;
   vector<float>   *tauEMRadius;
   vector<float>   *tauHadRadius;
   vector<float>   *tauEtEMAtEMScale;
   vector<float>   *tauEtHadAtEMScale;
   vector<float>   *tauIsolFrac;
   vector<float>   *tauCentFrac;
   vector<float>   *tauStripWidth2;
   vector<int>     *tauNumStrip;
   vector<float>   *tauTrkAvgDist;
   vector<float>   *tauDRmax;
   Float_t         metCaloPx;
   Float_t         metCaloPy;
   Float_t         metCaloET;
   Float_t         metCaloSumET;
   Float_t         metTrackPx;
   Float_t         metTrackPy;
   Float_t         metTrackET;
   Float_t         metTrackSumET;
   Float_t         metLocHadPx;
   Float_t         metLocHadPy;
   Float_t         metLocHadET;
   Float_t         metLocHadSumET;
   Float_t         metCoreLCTopoPx;
   Float_t         metCoreLCTopoPy;
   Float_t         metCoreLCTopoET;
   Float_t         metCoreLCTopoSumET;
   Float_t         metRefLCTopoPx;
   Float_t         metRefLCTopoPy;
   Float_t         metRefLCTopoET;
   Float_t         metRefLCTopoSumET;
   Float_t         metTruthPx;
   Float_t         metTruthPy;
   Float_t         metTruthET;
   Float_t         metTruthSumET;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_averagePU;   //!
   TBranch        *b_eventPU;   //!
   TBranch        *b_mcFlag;   //!
   TBranch        *b_numVtx;   //!
   TBranch        *b_pVtxX;   //!
   TBranch        *b_pVtxY;   //!
   TBranch        *b_pVtxZ;   //!
   TBranch        *b_pVtxXErr;   //!
   TBranch        *b_pVtxYErr;   //!
   TBranch        *b_pVtxZErr;   //!
   TBranch        *b_truthPVtxX;   //!
   TBranch        *b_truthPVtxY;   //!
   TBranch        *b_truthPVtxZ;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackTheta;   //!
   TBranch        *b_trackCharge;   //!
   TBranch        *b_trackqOverP;   //!
   TBranch        *b_trackD0;   //!
   TBranch        *b_trackZ0;   //!
   TBranch        *b_trackD0Err;   //!
   TBranch        *b_trackZ0Err;   //!
   TBranch        *b_trackqOverPErr;   //!
   TBranch        *b_trackDeltaZSinTheta;   //!
   TBranch        *b_trackMindR;   //!
   TBranch        *b_trackClass;   //!
   TBranch        *b_trackPassCut;   //!
   TBranch        *b_trackOutliers;   //!
   TBranch        *b_trackChiSqPerDof;   //!
   TBranch        *b_trackStdDevChi2OS;   //!
   TBranch        *b_trackNPixelHits;   //!
   TBranch        *b_trackNPixelHoles;   //!
   TBranch        *b_trackNPixelLayers;   //!
   TBranch        *b_trackNPixelOutliers;   //!
   TBranch        *b_nPixelGanged;   //!
   TBranch        *b_nPixelShared;   //!
   TBranch        *b_nPixelSplit;   //!
   TBranch        *b_nPixelDeadSensors;   //!
   TBranch        *b_trackNIBLUsedHits;   //!
   TBranch        *b_trackNBLUsedHits;   //!
   TBranch        *b_trackNSCTHits;   //!
   TBranch        *b_trackNSCTHoles;   //!
   TBranch        *b_nSCTOutliers;   //!
   TBranch        *b_nSCTDeadSensors;   //!
   TBranch        *b_trackNTRTHits;   //!
   TBranch        *b_nTRTOutliers;   //!
   TBranch        *b_nTRTHoles;   //!
   TBranch        *b_nTRTHTHits;   //!
   TBranch        *b_trackNSplitIBLHits;   //!
   TBranch        *b_trackNSplitBLHits;   //!
   TBranch        *b_trackNSharedIBLHits;   //!
   TBranch        *b_trackNSharedBLHits;   //!
   TBranch        *b_trackNIBLHits;   //!
   TBranch        *b_trackNBLHits;   //!
   TBranch        *b_trackNL1Hits;   //!
   TBranch        *b_trackNL2Hits;   //!
   TBranch        *b_trackNECHits;   //!
   TBranch        *b_trackNIBLExpectedHits;   //!
   TBranch        *b_trackNBLExpectedHits;   //!
   TBranch        *b_trackNSCTSharedHits;   //!
   TBranch        *b_trackTruncateddEdx;   //!
   TBranch        *b_trackNTruncateddEdx;   //!
   TBranch        *b_extrapolatedIBLX;   //!
   TBranch        *b_extrapolatedIBLY;   //!
   TBranch        *b_extrapolatedIBLZ;   //!
   TBranch        *b_extrapolatedBLX;   //!
   TBranch        *b_extrapolatedBLY;   //!
   TBranch        *b_extrapolatedBLZ;   //!
   TBranch        *b_extrapolatedL1X;   //!
   TBranch        *b_extrapolatedL1Y;   //!
   TBranch        *b_extrapolatedL1Z;   //!
   TBranch        *b_extrapolatedL2X;   //!
   TBranch        *b_extrapolatedL2Y;   //!
   TBranch        *b_extrapolatedL2Z;   //!
   TBranch        *b_truePt;   //!
   TBranch        *b_truePhi;   //!
   TBranch        *b_trueEta;   //!
   TBranch        *b_trueE;   //!
   TBranch        *b_trued0;   //!
   TBranch        *b_truez0;   //!
   TBranch        *b_truephi0;   //!
   TBranch        *b_truetheta;   //!
   TBranch        *b_trueqoverp;   //!
   TBranch        *b_truepdgid;   //!
   TBranch        *b_truebarcode;   //!
   TBranch        *b_truthmatchprob;   //!
   TBranch        *b_trueMindR;   //!
   TBranch        *b_genVtxR;   //!
   TBranch        *b_genVtxZ;   //!
   TBranch        *b_parentFlavour;   //!
   TBranch        *b_hitIsEndCap;   //!
   TBranch        *b_hitIsHole;   //!
   TBranch        *b_hitIsOutlier;   //!
   TBranch        *b_hitLayer;   //!
   TBranch        *b_hitEtaModule;   //!
   TBranch        *b_hitPhiModule;   //!
   TBranch        *b_hitCharge;   //!
   TBranch        *b_hitToT;   //!
   TBranch        *b_hitLVL1A;   //!
   TBranch        *b_hitNPixel;   //!
   TBranch        *b_hitNPixelX;   //!
   TBranch        *b_hitNPixelY;   //!
   TBranch        *b_hitBSerr;   //!
   TBranch        *b_hitDCSstate;   //!
   TBranch        *b_hitVBias;   //!
   TBranch        *b_hitTemp;   //!
   TBranch        *b_hitLorentzShift;   //!
   TBranch        *b_hitIsSplit;   //!
   TBranch        *b_hitIsEdge;   //!
   TBranch        *b_hitIsOverflow;   //!
   TBranch        *b_hitIsolation10x2;   //!
   TBranch        *b_hitIsolation20x4;   //!
   TBranch        *b_hitGlobalX;   //!
   TBranch        *b_hitGlobalY;   //!
   TBranch        *b_hitGlobalZ;   //!
   TBranch        *b_hitLocalX;   //!
   TBranch        *b_hitLocalY;   //!
   TBranch        *b_trkLocalX;   //!
   TBranch        *b_trkLocalY;   //!
   TBranch        *b_unbiasedResidualX;   //!
   TBranch        *b_unbiasedResidualY;   //!
   TBranch        *b_unbiasedPullX;   //!
   TBranch        *b_unbiasedPullY;   //!
   TBranch        *b_trkPhiOnSurface;   //!
   TBranch        *b_trkThetaOnSurface;   //!
   TBranch        *b_numTotalClustersPerModule;   //!
   TBranch        *b_numTotalPixelsPerModule;   //!
   TBranch        *b_RdoToT;   //!
   TBranch        *b_RdoCharge;   //!
   TBranch        *b_RdoPhi;   //!
   TBranch        *b_RdoEta;   //!
   TBranch        *b_hitTruthPdgId;   //!
   TBranch        *b_siHitPdgId;   //!
   TBranch        *b_siHitBarcode;   //!
   TBranch        *b_siHitStartPosX;   //!
   TBranch        *b_siHitStartPosY;   //!
   TBranch        *b_siHitEndPosX;   //!
   TBranch        *b_siHitEndPosY;   //!
   TBranch        *b_siHitEnergyDeposit;   //!
   TBranch        *b_truthEta;   //!
   TBranch        *b_truthPhi;   //!
   TBranch        *b_truthPt;   //!
   TBranch        *b_truthE;   //!
   TBranch        *b_truthCharge;   //!
   TBranch        *b_truthPdgId;   //!
   TBranch        *b_truthBarcode;   //!
   TBranch        *b_truthStatus;   //!
   TBranch        *b_muonEta;   //!
   TBranch        *b_muonPhi;   //!
   TBranch        *b_muonPt;   //!
   TBranch        *b_muonE;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPassJVT;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauM;   //!
   TBranch        *b_tauCharge;   //!
   TBranch        *b_tauNumTracksCharged;   //!
   TBranch        *b_tauNumTracksIsolation;   //!
   TBranch        *b_tauNumAllTracks;   //!
   TBranch        *b_tauTrackFilterProngs;   //!
   TBranch        *b_tauTrackFilterQuality;   //!
   TBranch        *b_tauNumClusters;   //!
   TBranch        *b_tauNumPi0s;   //!
   TBranch        *b_tauNumCells;   //!
   TBranch        *b_tauIsTauFlags;   //!
   TBranch        *b_tauBDTJetScore;   //!
   TBranch        *b_tauBDTEleScore;   //!
   TBranch        *b_tauEleMatchLikelihoodScore;   //!
   TBranch        *b_tauBDTJetScoreSigTrans;   //!
   TBranch        *b_tauRNNJetScore;   //!
   TBranch        *b_tauRNNJetScoreSigTrans;   //!
   TBranch        *b_tauIpZ0SinThetaSigLeadTrk;   //!
   TBranch        *b_tauEtOverPtLeadTrk;   //!
   TBranch        *b_tauLeadTrkPt;   //!
   TBranch        *b_tauIpSigLeadTrk;   //!
   TBranch        *b_tauMassTrkSys;   //!
   TBranch        *b_tauTrkWidth2;   //!
   TBranch        *b_tauTrFlightPathSig;   //!
   TBranch        *b_tauEMRadius;   //!
   TBranch        *b_tauHadRadius;   //!
   TBranch        *b_tauEtEMAtEMScale;   //!
   TBranch        *b_tauEtHadAtEMScale;   //!
   TBranch        *b_tauIsolFrac;   //!
   TBranch        *b_tauCentFrac;   //!
   TBranch        *b_tauStripWidth2;   //!
   TBranch        *b_tauNumStrip;   //!
   TBranch        *b_tauTrkAvgDist;   //!
   TBranch        *b_tauDRmax;   //!
   TBranch        *b_metCaloPx;   //!
   TBranch        *b_metCaloPy;   //!
   TBranch        *b_metCaloET;   //!
   TBranch        *b_metCaloSumET;   //!
   TBranch        *b_metTrackPx;   //!
   TBranch        *b_metTrackPy;   //!
   TBranch        *b_metTrackET;   //!
   TBranch        *b_metTrackSumET;   //!
   TBranch        *b_metLocHadPx;   //!
   TBranch        *b_metLocHadPy;   //!
   TBranch        *b_metLocHadET;   //!
   TBranch        *b_metLocHadSumET;   //!
   TBranch        *b_metCoreLCTopoPx;   //!
   TBranch        *b_metCoreLCTopoPy;   //!
   TBranch        *b_metCoreLCTopoET;   //!
   TBranch        *b_metCoreLCTopoSumET;   //!
   TBranch        *b_metRefLCTopoPx;   //!
   TBranch        *b_metRefLCTopoPy;   //!
   TBranch        *b_metRefLCTopoET;   //!
   TBranch        *b_metRefLCTopoSumET;   //!
   TBranch        *b_metTruthPx;   //!
   TBranch        *b_metTruthPy;   //!
   TBranch        *b_metTruthET;   //!
   TBranch        *b_metTruthSumET;   //!

   eventtree(TTree *tree=0);
   virtual ~eventtree();
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();

protected:
   virtual void initialize() {};
   virtual void execute() {};
   virtual void finalize() {};

};

#endif

#ifdef eventtree_cxx
eventtree::eventtree(TTree *tree) : fChain(0) 
{
   Init(tree);
}

eventtree::~eventtree()
{
}

void eventtree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trackPt = 0;
   trackPhi = 0;
   trackEta = 0;
   trackTheta = 0;
   trackCharge = 0;
   trackqOverP = 0;
   trackD0 = 0;
   trackZ0 = 0;
   trackD0Err = 0;
   trackZ0Err = 0;
   trackqOverPErr = 0;
   trackDeltaZSinTheta = 0;
   trackMindR = 0;
   trackClass = 0;
   trackPassCut = 0;
   trackOutliers = 0;
   trackChiSqPerDof = 0;
   trackStdDevChi2OS = 0;
   trackNPixelHits = 0;
   trackNPixelHoles = 0;
   trackNPixelLayers = 0;
   trackNPixelOutliers = 0;
   nPixelGanged = 0;
   nPixelShared = 0;
   nPixelSplit = 0;
   nPixelDeadSensors = 0;
   trackNIBLUsedHits = 0;
   trackNBLUsedHits = 0;
   trackNSCTHits = 0;
   trackNSCTHoles = 0;
   nSCTOutliers = 0;
   nSCTDeadSensors = 0;
   trackNTRTHits = 0;
   nTRTOutliers = 0;
   nTRTHoles = 0;
   nTRTHTHits = 0;
   trackNSplitIBLHits = 0;
   trackNSplitBLHits = 0;
   trackNSharedIBLHits = 0;
   trackNSharedBLHits = 0;
   trackNIBLHits = 0;
   trackNBLHits = 0;
   trackNL1Hits = 0;
   trackNL2Hits = 0;
   trackNECHits = 0;
   trackNIBLExpectedHits = 0;
   trackNBLExpectedHits = 0;
   trackNSCTSharedHits = 0;
   trackTruncateddEdx = 0;
   trackNTruncateddEdx = 0;
   extrapolatedIBLX = 0;
   extrapolatedIBLY = 0;
   extrapolatedIBLZ = 0;
   extrapolatedBLX = 0;
   extrapolatedBLY = 0;
   extrapolatedBLZ = 0;
   extrapolatedL1X = 0;
   extrapolatedL1Y = 0;
   extrapolatedL1Z = 0;
   extrapolatedL2X = 0;
   extrapolatedL2Y = 0;
   extrapolatedL2Z = 0;
   truePt = 0;
   truePhi = 0;
   trueEta = 0;
   trueE = 0;
   trued0 = 0;
   truez0 = 0;
   truephi0 = 0;
   truetheta = 0;
   trueqoverp = 0;
   truepdgid = 0;
   truebarcode = 0;
   truthmatchprob = 0;
   trueMindR = 0;
   genVtxR = 0;
   genVtxZ = 0;
   parentFlavour = 0;
   hitIsEndCap = 0;
   hitIsHole = 0;
   hitIsOutlier = 0;
   hitLayer = 0;
   hitEtaModule = 0;
   hitPhiModule = 0;
   hitCharge = 0;
   hitToT = 0;
   hitLVL1A = 0;
   hitNPixel = 0;
   hitNPixelX = 0;
   hitNPixelY = 0;
   hitBSerr = 0;
   hitDCSstate = 0;
   hitVBias = 0;
   hitTemp = 0;
   hitLorentzShift = 0;
   hitIsSplit = 0;
   hitIsEdge = 0;
   hitIsOverflow = 0;
   hitIsolation10x2 = 0;
   hitIsolation20x4 = 0;
   hitGlobalX = 0;
   hitGlobalY = 0;
   hitGlobalZ = 0;
   hitLocalX = 0;
   hitLocalY = 0;
   trkLocalX = 0;
   trkLocalY = 0;
   unbiasedResidualX = 0;
   unbiasedResidualY = 0;
   unbiasedPullX = 0;
   unbiasedPullY = 0;
   trkPhiOnSurface = 0;
   trkThetaOnSurface = 0;
   numTotalClustersPerModule = 0;
   numTotalPixelsPerModule = 0;
   RdoToT = 0;
   RdoCharge = 0;
   RdoPhi = 0;
   RdoEta = 0;
   hitTruthPdgId = 0;
   siHitPdgId = 0;
   siHitBarcode = 0;
   siHitStartPosX = 0;
   siHitStartPosY = 0;
   siHitEndPosX = 0;
   siHitEndPosY = 0;
   siHitEnergyDeposit = 0;
   truthEta = 0;
   truthPhi = 0;
   truthPt = 0;
   truthE = 0;
   truthCharge = 0;
   truthPdgId = 0;
   truthBarcode = 0;
   truthStatus = 0;
   muonEta = 0;
   muonPhi = 0;
   muonPt = 0;
   muonE = 0;
   jetEta = 0;
   jetPhi = 0;
   jetPt = 0;
   jetE = 0;
   jetPassJVT = 0;
   tauEta = 0;
   tauPhi = 0;
   tauPt = 0;
   tauM = 0;
   tauCharge = 0;
   tauNumTracksCharged = 0;
   tauNumTracksIsolation = 0;
   tauNumAllTracks = 0;
   tauTrackFilterProngs = 0;
   tauTrackFilterQuality = 0;
   tauNumClusters = 0;
   tauNumPi0s = 0;
   tauNumCells = 0;
   tauIsTauFlags = 0;
   tauBDTJetScore = 0;
   tauBDTEleScore = 0;
   tauEleMatchLikelihoodScore = 0;
   tauBDTJetScoreSigTrans = 0;
   tauRNNJetScore = 0;
   tauRNNJetScoreSigTrans = 0;
   tauIpZ0SinThetaSigLeadTrk = 0;
   tauEtOverPtLeadTrk = 0;
   tauLeadTrkPt = 0;
   tauIpSigLeadTrk = 0;
   tauMassTrkSys = 0;
   tauTrkWidth2 = 0;
   tauTrFlightPathSig = 0;
   tauEMRadius = 0;
   tauHadRadius = 0;
   tauEtEMAtEMScale = 0;
   tauEtHadAtEMScale = 0;
   tauIsolFrac = 0;
   tauCentFrac = 0;
   tauStripWidth2 = 0;
   tauNumStrip = 0;
   tauTrkAvgDist = 0;
   tauDRmax = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("averagePU", &averagePU, &b_averagePU);
   fChain->SetBranchAddress("eventPU", &eventPU, &b_eventPU);
   fChain->SetBranchAddress("mcFlag", &mcFlag, &b_mcFlag);
   fChain->SetBranchAddress("numVtx", &numVtx, &b_numVtx);
   fChain->SetBranchAddress("pVtxX", &pVtxX, &b_pVtxX);
   fChain->SetBranchAddress("pVtxY", &pVtxY, &b_pVtxY);
   fChain->SetBranchAddress("pVtxZ", &pVtxZ, &b_pVtxZ);
   fChain->SetBranchAddress("pVtxXErr", &pVtxXErr, &b_pVtxXErr);
   fChain->SetBranchAddress("pVtxYErr", &pVtxYErr, &b_pVtxYErr);
   fChain->SetBranchAddress("pVtxZErr", &pVtxZErr, &b_pVtxZErr);
   fChain->SetBranchAddress("truthPVtxX", &truthPVtxX, &b_truthPVtxX);
   fChain->SetBranchAddress("truthPVtxY", &truthPVtxY, &b_truthPVtxY);
   fChain->SetBranchAddress("truthPVtxZ", &truthPVtxZ, &b_truthPVtxZ);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackTheta", &trackTheta, &b_trackTheta);
   fChain->SetBranchAddress("trackCharge", &trackCharge, &b_trackCharge);
   fChain->SetBranchAddress("trackqOverP", &trackqOverP, &b_trackqOverP);
   fChain->SetBranchAddress("trackD0", &trackD0, &b_trackD0);
   fChain->SetBranchAddress("trackZ0", &trackZ0, &b_trackZ0);
   fChain->SetBranchAddress("trackD0Err", &trackD0Err, &b_trackD0Err);
   fChain->SetBranchAddress("trackZ0Err", &trackZ0Err, &b_trackZ0Err);
   fChain->SetBranchAddress("trackqOverPErr", &trackqOverPErr, &b_trackqOverPErr);
   fChain->SetBranchAddress("trackDeltaZSinTheta", &trackDeltaZSinTheta, &b_trackDeltaZSinTheta);
   fChain->SetBranchAddress("trackMindR", &trackMindR, &b_trackMindR);
   fChain->SetBranchAddress("trackClass", &trackClass, &b_trackClass);
   fChain->SetBranchAddress("trackPassCut", &trackPassCut, &b_trackPassCut);
   fChain->SetBranchAddress("trackOutliers", &trackOutliers, &b_trackOutliers);
   fChain->SetBranchAddress("trackChiSqPerDof", &trackChiSqPerDof, &b_trackChiSqPerDof);
   fChain->SetBranchAddress("trackStdDevChi2OS", &trackStdDevChi2OS, &b_trackStdDevChi2OS);
   fChain->SetBranchAddress("trackNPixelHits", &trackNPixelHits, &b_trackNPixelHits);
   fChain->SetBranchAddress("trackNPixelHoles", &trackNPixelHoles, &b_trackNPixelHoles);
   fChain->SetBranchAddress("trackNPixelLayers", &trackNPixelLayers, &b_trackNPixelLayers);
   fChain->SetBranchAddress("trackNPixelOutliers", &trackNPixelOutliers, &b_trackNPixelOutliers);
   fChain->SetBranchAddress("nPixelGanged", &nPixelGanged, &b_nPixelGanged);
   fChain->SetBranchAddress("nPixelShared", &nPixelShared, &b_nPixelShared);
   fChain->SetBranchAddress("nPixelSplit", &nPixelSplit, &b_nPixelSplit);
   fChain->SetBranchAddress("nPixelDeadSensors", &nPixelDeadSensors, &b_nPixelDeadSensors);
   fChain->SetBranchAddress("trackNIBLUsedHits", &trackNIBLUsedHits, &b_trackNIBLUsedHits);
   fChain->SetBranchAddress("trackNBLUsedHits", &trackNBLUsedHits, &b_trackNBLUsedHits);
   fChain->SetBranchAddress("trackNSCTHits", &trackNSCTHits, &b_trackNSCTHits);
   fChain->SetBranchAddress("trackNSCTHoles", &trackNSCTHoles, &b_trackNSCTHoles);
   fChain->SetBranchAddress("nSCTOutliers", &nSCTOutliers, &b_nSCTOutliers);
   fChain->SetBranchAddress("nSCTDeadSensors", &nSCTDeadSensors, &b_nSCTDeadSensors);
   fChain->SetBranchAddress("trackNTRTHits", &trackNTRTHits, &b_trackNTRTHits);
   fChain->SetBranchAddress("nTRTOutliers", &nTRTOutliers, &b_nTRTOutliers);
   fChain->SetBranchAddress("nTRTHoles", &nTRTHoles, &b_nTRTHoles);
   fChain->SetBranchAddress("nTRTHTHits", &nTRTHTHits, &b_nTRTHTHits);
   fChain->SetBranchAddress("trackNSplitIBLHits", &trackNSplitIBLHits, &b_trackNSplitIBLHits);
   fChain->SetBranchAddress("trackNSplitBLHits", &trackNSplitBLHits, &b_trackNSplitBLHits);
   fChain->SetBranchAddress("trackNSharedIBLHits", &trackNSharedIBLHits, &b_trackNSharedIBLHits);
   fChain->SetBranchAddress("trackNSharedBLHits", &trackNSharedBLHits, &b_trackNSharedBLHits);
   fChain->SetBranchAddress("trackNIBLHits", &trackNIBLHits, &b_trackNIBLHits);
   fChain->SetBranchAddress("trackNBLHits", &trackNBLHits, &b_trackNBLHits);
   fChain->SetBranchAddress("trackNL1Hits", &trackNL1Hits, &b_trackNL1Hits);
   fChain->SetBranchAddress("trackNL2Hits", &trackNL2Hits, &b_trackNL2Hits);
   fChain->SetBranchAddress("trackNECHits", &trackNECHits, &b_trackNECHits);
   fChain->SetBranchAddress("trackNIBLExpectedHits", &trackNIBLExpectedHits, &b_trackNIBLExpectedHits);
   fChain->SetBranchAddress("trackNBLExpectedHits", &trackNBLExpectedHits, &b_trackNBLExpectedHits);
   fChain->SetBranchAddress("trackNSCTSharedHits", &trackNSCTSharedHits, &b_trackNSCTSharedHits);
   fChain->SetBranchAddress("trackTruncateddEdx", &trackTruncateddEdx, &b_trackTruncateddEdx);
   fChain->SetBranchAddress("trackNTruncateddEdx", &trackNTruncateddEdx, &b_trackNTruncateddEdx);
   fChain->SetBranchAddress("extrapolatedIBLX", &extrapolatedIBLX, &b_extrapolatedIBLX);
   fChain->SetBranchAddress("extrapolatedIBLY", &extrapolatedIBLY, &b_extrapolatedIBLY);
   fChain->SetBranchAddress("extrapolatedIBLZ", &extrapolatedIBLZ, &b_extrapolatedIBLZ);
   fChain->SetBranchAddress("extrapolatedBLX", &extrapolatedBLX, &b_extrapolatedBLX);
   fChain->SetBranchAddress("extrapolatedBLY", &extrapolatedBLY, &b_extrapolatedBLY);
   fChain->SetBranchAddress("extrapolatedBLZ", &extrapolatedBLZ, &b_extrapolatedBLZ);
   fChain->SetBranchAddress("extrapolatedL1X", &extrapolatedL1X, &b_extrapolatedL1X);
   fChain->SetBranchAddress("extrapolatedL1Y", &extrapolatedL1Y, &b_extrapolatedL1Y);
   fChain->SetBranchAddress("extrapolatedL1Z", &extrapolatedL1Z, &b_extrapolatedL1Z);
   fChain->SetBranchAddress("extrapolatedL2X", &extrapolatedL2X, &b_extrapolatedL2X);
   fChain->SetBranchAddress("extrapolatedL2Y", &extrapolatedL2Y, &b_extrapolatedL2Y);
   fChain->SetBranchAddress("extrapolatedL2Z", &extrapolatedL2Z, &b_extrapolatedL2Z);
   fChain->SetBranchAddress("truePt", &truePt, &b_truePt);
   fChain->SetBranchAddress("truePhi", &truePhi, &b_truePhi);
   fChain->SetBranchAddress("trueEta", &trueEta, &b_trueEta);
   fChain->SetBranchAddress("trueE", &trueE, &b_trueE);
   fChain->SetBranchAddress("trued0", &trued0, &b_trued0);
   fChain->SetBranchAddress("truez0", &truez0, &b_truez0);
   fChain->SetBranchAddress("truephi0", &truephi0, &b_truephi0);
   fChain->SetBranchAddress("truetheta", &truetheta, &b_truetheta);
   fChain->SetBranchAddress("trueqoverp", &trueqoverp, &b_trueqoverp);
   fChain->SetBranchAddress("truepdgid", &truepdgid, &b_truepdgid);
   fChain->SetBranchAddress("truebarcode", &truebarcode, &b_truebarcode);
   fChain->SetBranchAddress("truthmatchprob", &truthmatchprob, &b_truthmatchprob);
   fChain->SetBranchAddress("trueMindR", &trueMindR, &b_trueMindR);
   fChain->SetBranchAddress("genVtxR", &genVtxR, &b_genVtxR);
   fChain->SetBranchAddress("genVtxZ", &genVtxZ, &b_genVtxZ);
   fChain->SetBranchAddress("parentFlavour", &parentFlavour, &b_parentFlavour);
   fChain->SetBranchAddress("hitIsEndCap", &hitIsEndCap, &b_hitIsEndCap);
   fChain->SetBranchAddress("hitIsHole", &hitIsHole, &b_hitIsHole);
   fChain->SetBranchAddress("hitIsOutlier", &hitIsOutlier, &b_hitIsOutlier);
   fChain->SetBranchAddress("hitLayer", &hitLayer, &b_hitLayer);
   fChain->SetBranchAddress("hitEtaModule", &hitEtaModule, &b_hitEtaModule);
   fChain->SetBranchAddress("hitPhiModule", &hitPhiModule, &b_hitPhiModule);
   fChain->SetBranchAddress("hitCharge", &hitCharge, &b_hitCharge);
   fChain->SetBranchAddress("hitToT", &hitToT, &b_hitToT);
   fChain->SetBranchAddress("hitLVL1A", &hitLVL1A, &b_hitLVL1A);
   fChain->SetBranchAddress("hitNPixel", &hitNPixel, &b_hitNPixel);
   fChain->SetBranchAddress("hitNPixelX", &hitNPixelX, &b_hitNPixelX);
   fChain->SetBranchAddress("hitNPixelY", &hitNPixelY, &b_hitNPixelY);
   fChain->SetBranchAddress("hitBSerr", &hitBSerr, &b_hitBSerr);
   fChain->SetBranchAddress("hitDCSstate", &hitDCSstate, &b_hitDCSstate);
   fChain->SetBranchAddress("hitVBias", &hitVBias, &b_hitVBias);
   fChain->SetBranchAddress("hitTemp", &hitTemp, &b_hitTemp);
   fChain->SetBranchAddress("hitLorentzShift", &hitLorentzShift, &b_hitLorentzShift);
   fChain->SetBranchAddress("hitIsSplit", &hitIsSplit, &b_hitIsSplit);
   fChain->SetBranchAddress("hitIsEdge", &hitIsEdge, &b_hitIsEdge);
   fChain->SetBranchAddress("hitIsOverflow", &hitIsOverflow, &b_hitIsOverflow);
   fChain->SetBranchAddress("hitIsolation10x2", &hitIsolation10x2, &b_hitIsolation10x2);
   fChain->SetBranchAddress("hitIsolation20x4", &hitIsolation20x4, &b_hitIsolation20x4);
   fChain->SetBranchAddress("hitGlobalX", &hitGlobalX, &b_hitGlobalX);
   fChain->SetBranchAddress("hitGlobalY", &hitGlobalY, &b_hitGlobalY);
   fChain->SetBranchAddress("hitGlobalZ", &hitGlobalZ, &b_hitGlobalZ);
   fChain->SetBranchAddress("hitLocalX", &hitLocalX, &b_hitLocalX);
   fChain->SetBranchAddress("hitLocalY", &hitLocalY, &b_hitLocalY);
   fChain->SetBranchAddress("trkLocalX", &trkLocalX, &b_trkLocalX);
   fChain->SetBranchAddress("trkLocalY", &trkLocalY, &b_trkLocalY);
   fChain->SetBranchAddress("unbiasedResidualX", &unbiasedResidualX, &b_unbiasedResidualX);
   fChain->SetBranchAddress("unbiasedResidualY", &unbiasedResidualY, &b_unbiasedResidualY);
   fChain->SetBranchAddress("unbiasedPullX", &unbiasedPullX, &b_unbiasedPullX);
   fChain->SetBranchAddress("unbiasedPullY", &unbiasedPullY, &b_unbiasedPullY);
   fChain->SetBranchAddress("trkPhiOnSurface", &trkPhiOnSurface, &b_trkPhiOnSurface);
   fChain->SetBranchAddress("trkThetaOnSurface", &trkThetaOnSurface, &b_trkThetaOnSurface);
   fChain->SetBranchAddress("numTotalClustersPerModule", &numTotalClustersPerModule, &b_numTotalClustersPerModule);
   fChain->SetBranchAddress("numTotalPixelsPerModule", &numTotalPixelsPerModule, &b_numTotalPixelsPerModule);
   fChain->SetBranchAddress("RdoToT", &RdoToT, &b_RdoToT);
   fChain->SetBranchAddress("RdoCharge", &RdoCharge, &b_RdoCharge);
   fChain->SetBranchAddress("RdoPhi", &RdoPhi, &b_RdoPhi);
   fChain->SetBranchAddress("RdoEta", &RdoEta, &b_RdoEta);
   fChain->SetBranchAddress("hitTruthPdgId", &hitTruthPdgId, &b_hitTruthPdgId);
   fChain->SetBranchAddress("siHitPdgId", &siHitPdgId, &b_siHitPdgId);
   fChain->SetBranchAddress("siHitBarcode", &siHitBarcode, &b_siHitBarcode);
   fChain->SetBranchAddress("siHitStartPosX", &siHitStartPosX, &b_siHitStartPosX);
   fChain->SetBranchAddress("siHitStartPosY", &siHitStartPosY, &b_siHitStartPosY);
   fChain->SetBranchAddress("siHitEndPosX", &siHitEndPosX, &b_siHitEndPosX);
   fChain->SetBranchAddress("siHitEndPosY", &siHitEndPosY, &b_siHitEndPosY);
   fChain->SetBranchAddress("siHitEnergyDeposit", &siHitEnergyDeposit, &b_siHitEnergyDeposit);
   fChain->SetBranchAddress("truthEta", &truthEta, &b_truthEta);
   fChain->SetBranchAddress("truthPhi", &truthPhi, &b_truthPhi);
   fChain->SetBranchAddress("truthPt", &truthPt, &b_truthPt);
   fChain->SetBranchAddress("truthE", &truthE, &b_truthE);
   fChain->SetBranchAddress("truthCharge", &truthCharge, &b_truthCharge);
   fChain->SetBranchAddress("truthPdgId", &truthPdgId, &b_truthPdgId);
   fChain->SetBranchAddress("truthBarcode", &truthBarcode, &b_truthBarcode);
   fChain->SetBranchAddress("truthStatus", &truthStatus, &b_truthStatus);
   fChain->SetBranchAddress("muonEta", &muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonPhi", &muonPhi, &b_muonPhi);
   fChain->SetBranchAddress("muonPt", &muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonE", &muonE, &b_muonE);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
   fChain->SetBranchAddress("jetPassJVT", &jetPassJVT, &b_jetPassJVT);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauM", &tauM, &b_tauM);
   fChain->SetBranchAddress("tauCharge", &tauCharge, &b_tauCharge);
   fChain->SetBranchAddress("tauNumTracksCharged", &tauNumTracksCharged, &b_tauNumTracksCharged);
   fChain->SetBranchAddress("tauNumTracksIsolation", &tauNumTracksIsolation, &b_tauNumTracksIsolation);
   fChain->SetBranchAddress("tauNumAllTracks", &tauNumAllTracks, &b_tauNumAllTracks);
   fChain->SetBranchAddress("tauTrackFilterProngs", &tauTrackFilterProngs, &b_tauTrackFilterProngs);
   fChain->SetBranchAddress("tauTrackFilterQuality", &tauTrackFilterQuality, &b_tauTrackFilterQuality);
   fChain->SetBranchAddress("tauNumClusters", &tauNumClusters, &b_tauNumClusters);
   fChain->SetBranchAddress("tauNumPi0s", &tauNumPi0s, &b_tauNumPi0s);
   fChain->SetBranchAddress("tauNumCells", &tauNumCells, &b_tauNumCells);
   fChain->SetBranchAddress("tauIsTauFlags", &tauIsTauFlags, &b_tauIsTauFlags);
   fChain->SetBranchAddress("tauBDTJetScore", &tauBDTJetScore, &b_tauBDTJetScore);
   fChain->SetBranchAddress("tauBDTEleScore", &tauBDTEleScore, &b_tauBDTEleScore);
   fChain->SetBranchAddress("tauEleMatchLikelihoodScore", &tauEleMatchLikelihoodScore, &b_tauEleMatchLikelihoodScore);
   fChain->SetBranchAddress("tauBDTJetScoreSigTrans", &tauBDTJetScoreSigTrans, &b_tauBDTJetScoreSigTrans);
   fChain->SetBranchAddress("tauRNNJetScore", &tauRNNJetScore, &b_tauRNNJetScore);
   fChain->SetBranchAddress("tauRNNJetScoreSigTrans", &tauRNNJetScoreSigTrans, &b_tauRNNJetScoreSigTrans);
   fChain->SetBranchAddress("tauIpZ0SinThetaSigLeadTrk", &tauIpZ0SinThetaSigLeadTrk, &b_tauIpZ0SinThetaSigLeadTrk);
   fChain->SetBranchAddress("tauEtOverPtLeadTrk", &tauEtOverPtLeadTrk, &b_tauEtOverPtLeadTrk);
   fChain->SetBranchAddress("tauLeadTrkPt", &tauLeadTrkPt, &b_tauLeadTrkPt);
   fChain->SetBranchAddress("tauIpSigLeadTrk", &tauIpSigLeadTrk, &b_tauIpSigLeadTrk);
   fChain->SetBranchAddress("tauMassTrkSys", &tauMassTrkSys, &b_tauMassTrkSys);
   fChain->SetBranchAddress("tauTrkWidth2", &tauTrkWidth2, &b_tauTrkWidth2);
   fChain->SetBranchAddress("tauTrFlightPathSig", &tauTrFlightPathSig, &b_tauTrFlightPathSig);
   fChain->SetBranchAddress("tauEMRadius", &tauEMRadius, &b_tauEMRadius);
   fChain->SetBranchAddress("tauHadRadius", &tauHadRadius, &b_tauHadRadius);
   fChain->SetBranchAddress("tauEtEMAtEMScale", &tauEtEMAtEMScale, &b_tauEtEMAtEMScale);
   fChain->SetBranchAddress("tauEtHadAtEMScale", &tauEtHadAtEMScale, &b_tauEtHadAtEMScale);
   fChain->SetBranchAddress("tauIsolFrac", &tauIsolFrac, &b_tauIsolFrac);
   fChain->SetBranchAddress("tauCentFrac", &tauCentFrac, &b_tauCentFrac);
   fChain->SetBranchAddress("tauStripWidth2", &tauStripWidth2, &b_tauStripWidth2);
   fChain->SetBranchAddress("tauNumStrip", &tauNumStrip, &b_tauNumStrip);
   fChain->SetBranchAddress("tauTrkAvgDist", &tauTrkAvgDist, &b_tauTrkAvgDist);
   fChain->SetBranchAddress("tauDRmax", &tauDRmax, &b_tauDRmax);
   fChain->SetBranchAddress("metCaloPx", &metCaloPx, &b_metCaloPx);
   fChain->SetBranchAddress("metCaloPy", &metCaloPy, &b_metCaloPy);
   fChain->SetBranchAddress("metCaloET", &metCaloET, &b_metCaloET);
   fChain->SetBranchAddress("metCaloSumET", &metCaloSumET, &b_metCaloSumET);
   fChain->SetBranchAddress("metTrackPx", &metTrackPx, &b_metTrackPx);
   fChain->SetBranchAddress("metTrackPy", &metTrackPy, &b_metTrackPy);
   fChain->SetBranchAddress("metTrackET", &metTrackET, &b_metTrackET);
   fChain->SetBranchAddress("metTrackSumET", &metTrackSumET, &b_metTrackSumET);
   fChain->SetBranchAddress("metLocHadPx", &metLocHadPx, &b_metLocHadPx);
   fChain->SetBranchAddress("metLocHadPy", &metLocHadPy, &b_metLocHadPy);
   fChain->SetBranchAddress("metLocHadET", &metLocHadET, &b_metLocHadET);
   fChain->SetBranchAddress("metLocHadSumET", &metLocHadSumET, &b_metLocHadSumET);
   fChain->SetBranchAddress("metCoreLCTopoPx", &metCoreLCTopoPx, &b_metCoreLCTopoPx);
   fChain->SetBranchAddress("metCoreLCTopoPy", &metCoreLCTopoPy, &b_metCoreLCTopoPy);
   fChain->SetBranchAddress("metCoreLCTopoET", &metCoreLCTopoET, &b_metCoreLCTopoET);
   fChain->SetBranchAddress("metCoreLCTopoSumET", &metCoreLCTopoSumET, &b_metCoreLCTopoSumET);
   fChain->SetBranchAddress("metRefLCTopoPx", &metRefLCTopoPx, &b_metRefLCTopoPx);
   fChain->SetBranchAddress("metRefLCTopoPy", &metRefLCTopoPy, &b_metRefLCTopoPy);
   fChain->SetBranchAddress("metRefLCTopoET", &metRefLCTopoET, &b_metRefLCTopoET);
   fChain->SetBranchAddress("metRefLCTopoSumET", &metRefLCTopoSumET, &b_metRefLCTopoSumET);
   fChain->SetBranchAddress("metTruthPx", &metTruthPx, &b_metTruthPx);
   fChain->SetBranchAddress("metTruthPy", &metTruthPy, &b_metTruthPy);
   fChain->SetBranchAddress("metTruthET", &metTruthET, &b_metTruthET);
   fChain->SetBranchAddress("metTruthSumET", &metTruthSumET, &b_metTruthSumET);
   Notify();
}

Bool_t eventtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef eventtree_cxx

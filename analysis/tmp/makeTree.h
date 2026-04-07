#include "eventtree.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

class makeTree : public eventtree {
  public:

    makeTree() {};
    makeTree(TTree *ntuple) : eventtree(ntuple) {};
    virtual ~makeTree() {};

    virtual void initialize(TString rout);
    virtual void execute();
    virtual void finalize();

  private:

    TFile * m_rootFile;
    TTree* m_trktree;

    int m_runNumber;
    int m_eventNumber;
    int m_lumiB;
    float m_averagePU;
    float m_eventPU;
    int m_mcFlag;
    float m_pileupWeight;
    float m_mcWeight;

    float m_pVtxX;
    float m_pVtxY;
    float m_pVtxZ;
    float m_pVtxXErr;
    float m_pVtxYErr;
    float m_pVtxZErr;
    float m_truthPVtxX;
    float m_truthPVtxY;
    float m_truthPVtxZ;

    int m_nPixelHits;
    int m_nPixelHoles;
    int m_nPixelLayers;
    int m_nPixelOutliers;
    int m_nIBLHits;
    int m_nSCTHits;
    int m_nSCTHoles;
    int m_nTRTHits;
    int m_nBLHits;
    int m_nIBLSplitHits;
    int m_nBLSplitHits;
    int m_nIBLSharedHits;
    int m_nBLSharedHits;
    int m_nL2Hits;
    int m_nL3Hits;
    int m_nECHits;

    int m_nIBLExpectedHits;
    int m_nBLExpectedHits;
    int m_nSCTSharedHits;

    float m_trackCharge;
    float m_trackPt;
    float m_trackPhi;
    float m_trackEta;
    float m_trackTheta;
    float m_trackqOverP;
    float m_trackD0;
    float m_trackZ0;
    float m_trackD0Err;
    float m_trackZ0Err;
    float m_trackqOverPErr;
    float m_trackDeltaZSinTheta;

    int m_trackClass;
    int m_trackPassCut;

    float m_pixeldEdx;
    int   m_nPixeldEdx;

    float m_trueD0;
    float m_trueZ0;
    float m_truePhi;
    float m_trueTheta;
    float m_trueqOverP;
    int m_truepdgid;
    int m_truebarcode;
    float m_truthMatchProb;

    float m_genVtxR;
    float m_genVtxZ;
    int m_parentFlav;

    float m_mindR;
    int m_hitSplitIBL;
    int m_hitSplitBL;
    int m_hitSplitL2;
    int m_hitSplitL3;

    std::vector<int> m_hitIsEndCap;
    std::vector<int> m_hitIsHole;
    std::vector<int> m_hitIsOutlier;
    std::vector<int> m_hitIs3D;
    std::vector<int> m_hitIsIBL;
    std::vector<int> m_hitLayer; 
    std::vector<float> m_hitCharge; 
    std::vector<int> m_hitToT; 
    std::vector<int> m_hitLVL1A; 
    std::vector<int> m_hitNPixelX; 
    std::vector<int> m_hitNPixelY; 
    std::vector<int> m_hitEtaModule;  
    std::vector<int> m_hitPhiModule;
    std::vector<float> m_hitVBias; 
    std::vector<float> m_hitVDep; 
    std::vector<float> m_hitTemp; 
    std::vector<float> m_hitLorentzShift; 
    std::vector<int> m_hitIsSplit;
    std::vector<float> m_hitGlobalX;   
    std::vector<float> m_hitGlobalY;   
    std::vector<float> m_hitGlobalZ;   
    std::vector<float> m_hitLocalX; 
    std::vector<float> m_hitLocalY; 

    std::vector<float> m_g4LocalX; 
    std::vector<float> m_g4LocalY; 
    std::vector<int>   m_g4Barcode; 
    std::vector<int>   m_g4PdgId; 
    std::vector<float> m_g4EnergyDeposit; 

    std::vector<float> m_trkLocalX; 
    std::vector<float> m_trkLocalY; 
    std::vector<float> m_hitLocalErrorX; 
    std::vector<float> m_hitLocalErrorY; 
    std::vector<float> m_trkLocalErrorX; 
    std::vector<float> m_trkLocalErrorY; 
    std::vector<float> m_unbiasedResidualX; 
    std::vector<float> m_unbiasedResidualY; 
    std::vector<float> m_biasedResidualX; 
    std::vector<float> m_biasedResidualY;
    std::vector<float> m_unbiasedPullX; 
    std::vector<float> m_unbiasedPullY; 
    std::vector<float> m_biasedPullX; 
    std::vector<float> m_biasedPullY;

    std::vector<float> m_minTrkG4DistX; 
    std::vector<float> m_minTrkG4DistY; 

    std::vector<float> m_minTrkDistX; 
    std::vector<float> m_minTrkDistY; 

    std::vector<float> m_trkPhiOnSurface; 
    std::vector<float> m_trkThetaOnSurface; 

    std::vector< int > m_clus_nContributingPtcs;
    std::vector< int > m_clus_nContributingPU;
    std::vector< int > m_clus_contributingPtcsBarcode;

    std::vector< int > m_clus_pixelIBLIndex;
    std::vector< int > m_clus_pixelIBLRow;
    std::vector< int > m_clus_pixelIBLCol;
    std::vector< float > m_clus_pixelIBLCharge;
    std::vector< int > m_clus_pixelIBLToT;

    std::vector< int > m_clus_pixelBLIndex;
    std::vector< int > m_clus_pixelBLRow;
    std::vector< int > m_clus_pixelBLCol;
    std::vector< float > m_clus_pixelBLCharge;
    std::vector< int > m_clus_pixelBLToT;

    std::vector< int > m_clus_pixelL2Index;
    std::vector< int > m_clus_pixelL2Row;
    std::vector< int > m_clus_pixelL2Col;
    std::vector< float > m_clus_pixelL2Charge;
    std::vector< int > m_clus_pixelL2ToT;

    std::vector< int > m_clus_pixelL3Index;
    std::vector< int > m_clus_pixelL3Row;
    std::vector< int > m_clus_pixelL3Col;
    std::vector< float > m_clus_pixelL3Charge;
    std::vector< int > m_clus_pixelL3ToT;

};


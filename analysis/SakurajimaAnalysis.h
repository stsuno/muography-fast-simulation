#include "track_tree.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"

class SakurajimaAnalysis : public track_tree {
  public:

    SakurajimaAnalysis() {};
    SakurajimaAnalysis(TTree *ntuple) : track_tree(ntuple) {};
    virtual ~SakurajimaAnalysis() {};

    virtual void initialize(TString rout);
    virtual void execute();
    virtual void finalize();

  private:
    TH1F* hist_track_nhits;
    TH1F* hist_track_chi2;
    TH1F* hist_track_C;
    TH1F* hist_track_D;
    TH1F* hist_track_phi;
    TH1F* hist_track_theta;

    TH2F* hist2D_mountain;
    TH2F* hist2D_sky;

    TFile * m_rootFile;
};


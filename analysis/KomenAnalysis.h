#include "track_tree.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"

class KomenAnalysis : public track_tree {
  public:

    KomenAnalysis() {};
    KomenAnalysis(TTree *ntuple) : track_tree(ntuple) {};
    virtual ~KomenAnalysis() {};

    virtual void initialize(TString rout);
    virtual void execute();
    virtual void finalize();

  private:

    TH1F* hist_track_theta;
    TH1F* hist_track_phi;

    TFile * m_rootFile;
};


#include "track_tree.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"

class DemoTest : public track_tree {
  public:

    DemoTest() {};
    DemoTest(TTree *ntuple) : track_tree(ntuple) {};
    virtual ~DemoTest() {};

    virtual void initialize(TString rout);
    virtual void execute();
    virtual void finalize();

  private:
    TH1F* hist_track_nhits;
    TH1F* hist_track_chi2;
    TH1F* hist_track_C;
    TH1F* hist_track_D;
    TH1F* hist_track_errA;
    TH1F* hist_track_errB;
    TH1F* hist_track_errC;
    TH1F* hist_track_errD;
    TH1F* hist_track_phi;
    TH1F* hist_track_theta;

    TH1F* hist_muon_energy;
    TH1F* hist_muon_loss;
    TH1F* hist_muon_detenergy;
    TH1F* hist_muon_charge;
    TH1F* hist_muon_phi;
    TH1F* hist_muon_theta;
    TH1F* hist_diff_phi;
    TH1F* hist_diff_theta;

    TH1F* hist_residual_x1;
    TH1F* hist_residual_x2;
    TH1F* hist_residual_x3;
    TH1F* hist_residual_x4;

    TH1F* hist_residual_y1;
    TH1F* hist_residual_y2;
    TH1F* hist_residual_y3;
    TH1F* hist_residual_y4;

    TH2F* hist2D_xy50;
    TH2F* hist2D_xy100;

    TFile * m_rootFile;
};


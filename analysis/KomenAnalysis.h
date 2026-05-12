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

    void SetMode(int mode) { m_mode=mode; }
    void SetCheckFile(std::string file) { m_checkfile=file; }
    void SetDirectryName(std::string directry) { m_FDirectoryName=directry; }
    void SetMWPCPosition(TVector3 position) { m_positionMWPC=position; }

    virtual void initialize(TString rout);
    virtual void execute();
    virtual void finalize();

  private:
    int m_mode;
    TVector3 m_positionMWPC;
    std::string m_checkfile;
    std::string m_FDirectoryName;

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
    TH1F* hist_narrow_theta;

    TH1F* hist_mb_phi;
    TH1F* hist_mb_theta;
    TH1F* diff_mb_phi;
    TH1F* diff_mb_theta;
    TH1F* diff_mb_deltaR;
    TH1F* diff_mb_C;
    TH1F* diff_mb_D;
    TH1F* diff_mb_dist;

    TH1F* hist_mc_phi;
    TH1F* hist_mc_theta;
    TH1F* diff_mc_phi;
    TH1F* diff_mc_theta;
    TH1F* diff_mc_deltaR;
    TH1F* diff_mc_C;
    TH1F* diff_mc_D;
    TH1F* diff_mc_dist;

    TH1F* hist_ma_phi;
    TH1F* hist_ma_theta;
    TH1F* diff_ma_phi;
    TH1F* diff_ma_theta;
    TH1F* diff_ma_deltaR;
    TH1F* diff_ma_C;
    TH1F* diff_ma_D;
    TH1F* diff_ma_dist;

    TH1F* hist_bc_phi;
    TH1F* hist_bc_theta;

    TH1F* hist_mb_simpletheta;
    TH1F* hist_mc_simpletheta;
    TH1F* hist_ma_simpletheta;
    TH1F* hist_bc_simpletheta;

    TFile * m_rootFile;
};


#ifdef __CLING__
R__ADD_INCLUDE_PATH(.);
#endif

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TEnv.h>
#include <TMath.h>
#include <TError.h>
#include <TSystem.h>

#include "Simulation/SciBar100x100.h"

#include "Reconstruction/ClusterMaker.h"
#include "Reconstruction/LineFitter.h"
#include "Reconstruction/TrackFinder.h"
#include "Reconstruction/Track.h"

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

int ReconstructionTest(bool withAnomaly=true) {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

  std::string output_dir = "data/output/test";

//  bool withAnomaly = false;
  std::string datafile = "test_anomaly.dat";
  std::string rootfile = "test_anomaly.root";
  std::string pdffile  = "recotest_anomaly.pdf";

  if (!withAnomaly) {
    datafile = "test_noanomaly.dat";
    rootfile = "test_noanomaly.root";
    pdffile  = "recotest_noanomaly.pdf";
  }

  // 2. Setup Output
  std::string output_path = output_dir + "/" + rootfile;
  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    Info("Reconstruction", "Created output directory: %s", output_dir.c_str());
  }
  TFile* output_root = new TFile(output_path.c_str(), "RECREATE");

  // 3. Define TTree
  TTree* track_tree = new TTree("track_tree", "Reconstructed Tracks");
  double m_A, m_B, m_C, m_D;
  double m_Aerr, m_Berr, m_Cerr, m_Derr;
  double m_chisq, m_p_value, m_phi, m_theta;
  int m_num_hits;
  double m_muene, m_muloss, m_muphi, m_mutheta;
  int m_muchg;
  std::vector<double>* m_residual_x = new std::vector<double>();
  std::vector<double>* m_residual_y = new std::vector<double>();
  track_tree->Branch("num_hits",      &m_num_hits,      "num_hits/I");
  track_tree->Branch("chisq",         &m_chisq,         "chisq/D");
  track_tree->Branch("p_value",       &m_p_value,       "p_value/D");
  track_tree->Branch("A",             &m_A,             "A/D");  // Slope dX/dZ
  track_tree->Branch("B",             &m_B,             "B/D");  // Slope dY/dZ
  track_tree->Branch("C",             &m_C,             "C/D");  // Offset X at Z=0
  track_tree->Branch("D",             &m_D,             "D/D");  // Offset Y at Z=0
  track_tree->Branch("errorA",        &m_Aerr,          "errorA/D");
  track_tree->Branch("errorB",        &m_Berr,          "errorB/D");
  track_tree->Branch("errorC",        &m_Cerr,          "errorC/D");
  track_tree->Branch("errorD",        &m_Derr,          "errorD/D");
  track_tree->Branch("phi",           &m_phi,           "phi/D");
  track_tree->Branch("theta",         &m_theta,         "theta/D");
  track_tree->Branch("residual_x", "std::vector<double>", &m_residual_x);
  track_tree->Branch("residual_y", "std::vector<double>", &m_residual_y);
  track_tree->Branch("muon_energy",   &m_muene,         "muon_energy/D");
  track_tree->Branch("muon_loss",     &m_muloss,        "muon_loss/D");
  track_tree->Branch("muon_charge",   &m_muchg,         "muon_charge/I");
  track_tree->Branch("muon_phi",      &m_muphi,         "muon_phi/D");
  track_tree->Branch("muon_theta",    &m_mutheta,       "muon_theta/D");

  // 4. Open Input Data
  std::string input_path = output_dir + "/" + datafile;
  std::ifstream input_dat(input_path.c_str());
  if (!input_dat.is_open()) {
    SysError("Reconstruction", "Failed to open input file: %s", input_path.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  // 5. Initialize Reconstruction Tools
  SciBar100x100* detector = new SciBar100x100(0.0,0.0,0.0);
  detector -> SetLineColor(3);
  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  ClusterMaker<SciBar100x100*> cluster_maker(detector);
  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::CombinedLocal);
//  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::Global);

  TrackFinder track_finder(detector->GetNumLayers(), 1000000.0);
  track_finder.SetMinimumHits(3);

//  TrackFinder track_finder(8, 1000000.0);
//  track_finder.SetMinimumHits(6);

  std::map<std::string, double> start_values;
  start_values["A"] = 0.0;
  start_values["B"] = 0.0;
  start_values["C"] = 0.0;
  start_values["D"] = 0.0;
  track_finder.SetStartValues(start_values);

  // 6. Event Loop
  int n_events = 0;
  int i_events = 0;
  std::string line;

  // Preliminary count for progress bar
  while (std::getline(input_dat, line)) n_events++;
  input_dat.clear();
  input_dat.seekg(0, std::ios::beg);

  if (n_events == 0) {
    Warning("Reconstruction", "Input file is empty: %s", input_path.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  int n = std::to_string(n_events).length();
  int step_events = static_cast<int>(TMath::Power(10, TMath::Max(0, n-1)));

  while (std::getline(input_dat, line)) {
    if (i_events % step_events == 0) { Info("Reconstruction", "Progress: %*d / %d processed", n, i_events, n_events); }
    i_events++;

    if (line.empty()) continue;

    std::stringstream data(line);
    double muene, muloss, muphi, mutheta;
    int muchg;
    data >> muene >> muloss >> muchg >> mutheta >> muphi;
    int hit_id;
    std::vector<int> hit_ids;
    while (data >> hit_id) hit_ids.push_back(hit_id);

    // Step A: Clustering (Raw hits -> SpacePoints)
    std::vector<SpacePoint> all_space_points = cluster_maker.Execute(hit_ids);
    if (g_debug_level>0) {
     Info("Reconstruction", "Event %d: SpacePoints = %lu", i_events, all_space_points.size());
    }

    // Step B: Track Finding (Combinatorial Search and Fitting)
    track_finder.GetLineFit().SetFittingParameters(&start_values, nullptr, nullptr);
    Track best_track = track_finder.GetBestTrack(all_space_points);

    // Step C: Recording (Fill TTree if a valid track was found)
    if (best_track.IsValid()) {
      m_num_hits      = best_track.GetNumHits();
      m_chisq         = best_track.GetChiSquare();
      m_p_value       = best_track.GetPValue();
      m_A             = best_track.GetA();
      m_B             = best_track.GetB();
      m_C             = best_track.GetC();
      m_D             = best_track.GetD();
      m_Aerr          = best_track.GetErrorA();
      m_Berr          = best_track.GetErrorB();
      m_Cerr          = best_track.GetErrorC();
      m_Derr          = best_track.GetErrorD();
      m_phi           = best_track.GetPhi()*TMath::RadToDeg();
      m_theta         = best_track.GetTheta()*TMath::RadToDeg();
      *m_residual_x   = best_track.GetResidualX();
      *m_residual_y   = best_track.GetResidualY();
      m_muene         = muene;
      m_muloss        = muloss;
      m_muchg         = muchg;
      m_muphi         = muphi*TMath::RadToDeg()>180.0?muphi*TMath::RadToDeg()-360.0:muphi*TMath::RadToDeg();
      m_mutheta       = mutheta*TMath::RadToDeg();
      track_tree->Fill();
    }

  }

  Info("Reconstruction", "Progress: %*d / %d completed", n, n_events, n_events);
  output_root->Write();
  output_root->Close();
  Info("Reconstruction", "Saved histograms to %s", output_path.c_str());

  return 0;
}


#ifndef __CLING__
int main(int argc, char** argv) {
  std::string chk = argv[1];
  bool withAnomaly=true;
  if (chk.find("false")!=std::string::npos) { withAnomaly=false; }
  return ReconstructionTest(withAnomaly);
}
#endif

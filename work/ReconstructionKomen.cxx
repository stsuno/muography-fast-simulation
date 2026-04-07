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

int ReconstructionKomen(std::string model="ALL") {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

  std::string output_dir = "data/output/komen";

  std::string datafileA    = "komen_singleA.dat";
  std::string datafileB    = "komen_singleB.dat";
  std::string datafileM    = "komen_singleM.dat";
  std::string datafileMB_M = "komen_doubleMB_M.dat";
  std::string datafileMB_B = "komen_doubleMB_B.dat";
  std::string datafileMA_M = "komen_doubleMA_M.dat";
  std::string datafileMA_A = "komen_doubleMA_A.dat";

  std::string rootfileA    = "komen_singleA.root";
  std::string rootfileB    = "komen_singleB.root";
  std::string rootfileM    = "komen_singleM.root";
  std::string rootfileMB_M = "komen_doubleMB_M.root";
  std::string rootfileMB_B = "komen_doubleMB_B.root";
  std::string rootfileMA_M = "komen_doubleMA_M.root";
  std::string rootfileMA_A = "komen_doubleMA_A.root";

  std::string pdffile      = "komen.pdf";

  TVector3 detA_position(0.0,0.0,0.0);
  TVector3 detB_position(2000.0,0.0,0.0);
  TVector3 MWPC_position(2000.0,0.0,2000.0);

  // 2. Setup Output
//  std::string output_path = output_dir + "/" + rootfile;
  std::string output_path = output_dir + "/";
  if (model=="SingleA")    { output_path+=rootfileA; }
  if (model=="SingleB")    { output_path+=rootfileB; }
  if (model=="SingleMWPC") { output_path+=rootfileM; }
  if (model=="DoubleMB_M") { output_path+=rootfileMB_M; }
  if (model=="DoubleMB_B") { output_path+=rootfileMB_B; }
  if (model=="DoubleMA_M") { output_path+=rootfileMA_M; }
  if (model=="DoubleMA_A") { output_path+=rootfileMA_A; }

  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    Info("Reconstruction", "Created output directory: %s", output_dir.c_str());
  }
  TFile* output_root = new TFile(output_path.c_str(), "RECREATE");

  // 3. Define TTree
  TTree* track_tree = new TTree("track_tree", "Reconstructed Tracks");
  double m_A, m_B, m_C, m_D;
  double m_Aerr, m_Berr, m_Cerr, m_Derr;
  double m_chisq, m_reduced_chisq, m_p_value, m_phi, m_theta;
  int m_ndf, m_num_hits;
  std::vector<double>* m_residual_x = new std::vector<double>();
  std::vector<double>* m_residual_y = new std::vector<double>();
  track_tree->Branch("num_hits",      &m_num_hits,      "num_hits/I");
  track_tree->Branch("chisq",         &m_chisq,         "chisq/D");
  track_tree->Branch("ndf",           &m_ndf,           "ndf/I");
  track_tree->Branch("reduced_chisq", &m_reduced_chisq, "reduced_chisq/D");
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

  // 4. Open Input Data
  std::string input_path = output_dir + "/";
//  std::string input_path = output_dir + "/" + datafile;

  if (model=="SingleA")    { input_path+=datafileA; }
  if (model=="SingleB")    { input_path+=datafileB; }
  if (model=="SingleMWPC") { input_path+=datafileM; }
  if (model=="DoubleMB_M") { input_path+=datafileMB_M; }
  if (model=="DoubleMB_B") { input_path+=datafileMB_B; }
  if (model=="DoubleMA_M") { input_path+=datafileMA_M; }
  if (model=="DoubleMA_A") { input_path+=datafileMA_A; }

  std::ifstream input_dat(input_path.c_str());
  if (!input_dat.is_open()) {
    SysError("Reconstruction", "Failed to open input file: %s", input_path.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  // 5. Initialize Reconstruction Tools
  SciBar100x100* detector = new SciBar100x100(0.0,0.0,0.0);
  TVector3 detector_position;

  if (model=="SingleA")    { detector_position = detA_position; }
  if (model=="SingleB")    { detector_position = detB_position; }
  if (model=="SingleMWPC") { detector_position = MWPC_position; }
  if (model=="DoubleMB_M") { detector_position = MWPC_position; }
  if (model=="DoubleMB_B") { detector_position = detB_position; }
  if (model=="DoubleMA_M") { detector_position = MWPC_position; }
  if (model=="DoubleMA_A") { detector_position = detA_position; }

  detector -> MoveXYZ(detector_position.X(),detector_position.Y(),detector_position.Z());
  detector -> SetLineColor(3);
  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  ClusterMaker<SciBar100x100*> cluster_maker(detector);
  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::Global);
//  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::CombinedLocal);
//  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::CombinedGlobal);

  TrackFinder track_finder(8, 1000000.0);
  track_finder.SetMinimumHits(6);

//  TrackFinder track_finder(4, 1000000.0);
//  track_finder.SetMinimumHits(3);
//  track_finder.SetMinimumHits(4);

//  std::map<std::string, double> start_values;
//  start_values["A"] = 0.0;
//  start_values["B"] = 0.0;
//  start_values["C"] = detector_position.X();
//  start_values["D"] = detector_position.Y();
//  track_finder.SetStartValues(start_values);

  // 6. Event Loop
  int n_events = 0;
  int i_events = 0;
  std::string line;

  // Preliminary count for progress bar

  while (std::getline(input_dat, line)) {
    n_events++;
  }
  input_dat.clear();
  input_dat.seekg(0, std::ios::beg);

  if (n_events == 0) {
    Warning("Reconstruction", "Input file is empty: %s", input_path.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  int n = std::to_string(n_events).length();
  int step_events = static_cast<int>(TMath::Power(10, TMath::Max(0, n-1)));

  while (std::getline(input_dat, line)) {
    if (i_events%step_events==0) { 
      Info("Reconstruction", "Progress: %*d / %d processed", n, i_events, n_events); 
    }
    i_events++;

    if (line.empty()) continue;

    std::stringstream data(line);
    int hit_id;
    std::vector<int> hit_ids;
    while (data >> hit_id) { hit_ids.push_back(hit_id); } 

    // Step A: Clustering (Raw hits -> SpacePoints)
    std::vector<SpacePoint> all_space_points = cluster_maker.Execute(hit_ids);
    if (g_debug_level>0) {
     Info("Reconstruction", "Event %d: SpacePoints = %lu", i_events, all_space_points.size());
    }

    // Step A': Shift the space point to the reference point
    std::vector<SpacePoint> space_points;
    for (const auto& space_point : all_space_points) {
      SpacePoint sp(space_point.X()-detector_position.X(),
                    space_point.Y()-detector_position.Y(),
                    space_point.Z()-detector_position.Z());
      sp.SetErrorX(space_point.GetErrorX());
      sp.SetErrorY(space_point.GetErrorY());
      sp.SetErrorZ(space_point.GetErrorZ());
      sp.SetUnitIndex(space_point.GetUnitIndex());

      space_points.push_back(sp);

//      std::cout << sp.X() << " " << sp.Y() << " " << sp.Z() << " " 
//                << sp.GetErrorX() << " " << sp.GetErrorY() << " " << sp.GetErrorZ() << " " 
//                << space_point.GetUnitIndex() << std::endl;


    }

    // Step B: Track Finding (Combinatorial Search and Fitting)
    Track best_track = track_finder.GetBestTrack(space_points);
//    Track best_track = track_finder.GetBestTrack(all_space_points);

    // Step C: Recording (Fill TTree if a valid track was found)
    if (best_track.IsValid()) {

      m_num_hits      = best_track.GetNumHits();
      m_chisq         = best_track.GetChiSquare();
      m_ndf           = best_track.GetNDF();
      m_reduced_chisq = best_track.GetReducedChiSquare();
      m_p_value       = best_track.GetPValue();
      m_A             = best_track.GetA();
      m_B             = best_track.GetB();
      m_C             = best_track.GetC();
      m_D             = best_track.GetD();
      m_Aerr          = best_track.GetErrorA();
      m_Berr          = best_track.GetErrorB();
      m_Cerr          = best_track.GetErrorC();
      m_Derr          = best_track.GetErrorD();
      m_phi           = best_track.GetPhi() * TMath::RadToDeg();
      m_theta         = best_track.GetTheta() * TMath::RadToDeg();
      *m_residual_x   = best_track.GetResidualX();
      *m_residual_y   = best_track.GetResidualY();

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
  std::string model = argv[1];
  return ReconstructionKomen(model);
}
#endif

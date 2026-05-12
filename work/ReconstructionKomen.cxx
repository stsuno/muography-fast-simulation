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

  std::string output_dir = "data/output/komen_g00pm20p";
//  std::string output_dir = "data/output/komen_g00pm500p";
//  std::string output_dir = "data/output/komen_nominal";

  std::string pdffile      = "komen.pdf";

  TVector3 detA_position(0.0,0.0,0.0);
  TVector3 detB_position(2000.0,0.0,0.0);
  TVector3 detC_position(1000.0,0.0,200.0);
  TVector3 MWPC_position(2000.0,0.0,2400.0);

  // 2. Setup Output
//  std::string output_path = output_dir + "/" + rootfile;
  std::string output_path = output_dir + "/";
  if (model=="SingleMWPC")  { output_path+="komen_singleM.root"; }
  if (model=="SingleB")     { output_path+="komen_singleB.root"; }
  if (model=="SingleC")     { output_path+="komen_singleC.root"; }
  if (model=="SingleA")     { output_path+="komen_singleA.root"; }
  if (model=="DoubleM1B_M") { output_path+="komen_doubleM1B_M.root"; }
  if (model=="DoubleM1B_B") { output_path+="komen_doubleM1B_B.root"; }
  if (model=="DoubleM2B_M") { output_path+="komen_doubleM2B_M.root"; }
  if (model=="DoubleM2B_B") { output_path+="komen_doubleM2B_B.root"; }
  if (model=="DoubleM3B_M") { output_path+="komen_doubleM3B_M.root"; }
  if (model=="DoubleM3B_B") { output_path+="komen_doubleM3B_B.root"; }
  if (model=="DoubleM4B_M") { output_path+="komen_doubleM4B_M.root"; }
  if (model=="DoubleM4B_B") { output_path+="komen_doubleM4B_B.root"; }
  if (model=="DoubleM5B_M") { output_path+="komen_doubleM5B_M.root"; }
  if (model=="DoubleM5B_B") { output_path+="komen_doubleM5B_B.root"; }
  if (model=="DoubleM1C_M") { output_path+="komen_doubleM1C_M.root"; }
  if (model=="DoubleM1C_C") { output_path+="komen_doubleM1C_C.root"; }
  if (model=="DoubleM2C_M") { output_path+="komen_doubleM2C_M.root"; }
  if (model=="DoubleM2C_C") { output_path+="komen_doubleM2C_C.root"; }
  if (model=="DoubleM3C_M") { output_path+="komen_doubleM3C_M.root"; }
  if (model=="DoubleM3C_C") { output_path+="komen_doubleM3C_C.root"; }
  if (model=="DoubleM4C_M") { output_path+="komen_doubleM4C_M.root"; }
  if (model=="DoubleM4C_C") { output_path+="komen_doubleM4C_C.root"; }
  if (model=="DoubleM5C_M") { output_path+="komen_doubleM5C_M.root"; }
  if (model=="DoubleM5C_C") { output_path+="komen_doubleM5C_C.root"; }
  if (model=="DoubleM1A_M") { output_path+="komen_doubleM1A_M.root"; }
  if (model=="DoubleM1A_A") { output_path+="komen_doubleM1A_A.root"; }
  if (model=="DoubleM2A_M") { output_path+="komen_doubleM2A_M.root"; }
  if (model=="DoubleM2A_A") { output_path+="komen_doubleM2A_A.root"; }
  if (model=="DoubleM3A_M") { output_path+="komen_doubleM3A_M.root"; }
  if (model=="DoubleM3A_A") { output_path+="komen_doubleM3A_A.root"; }
  if (model=="DoubleM4A_M") { output_path+="komen_doubleM4A_M.root"; }
  if (model=="DoubleM4A_A") { output_path+="komen_doubleM4A_A.root"; }
  if (model=="DoubleM5A_M") { output_path+="komen_doubleM5A_M.root"; }
  if (model=="DoubleM5A_A") { output_path+="komen_doubleM5A_A.root"; }
  if (model=="DoubleBC_B")  { output_path+="komen_doubleBC_B.root"; }
  if (model=="DoubleBC_C")  { output_path+="komen_doubleBC_C.root"; }
  if (model=="DoubleCA_C")  { output_path+="komen_doubleCA_C.root"; }
  if (model=="DoubleCA_A")  { output_path+="komen_doubleCA_A.root"; }
  if (model=="DoubleBC_BC") { output_path+="komen_doubleBC_BC.root"; }

  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    Info("Reconstruction", "Created output directory: %s", output_dir.c_str());
  }
  TFile* output_root = new TFile(output_path.c_str(), "RECREATE");

  // 3. Define TTree
  TTree* track_tree = new TTree("track_tree", "Reconstructed Tracks");
  int m_counter;
  double m_A, m_B, m_C, m_D;
  double m_Aerr, m_Berr, m_Cerr, m_Derr;
  double m_chisq, m_reduced_chisq, m_p_value, m_phi, m_theta;
  int m_ndf, m_num_hits;
  std::vector<double>* m_residual_x = new std::vector<double>();
  std::vector<double>* m_residual_y = new std::vector<double>();
  track_tree->Branch("event_number",  &m_counter,       "event_number/I");
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

  if (model=="SingleMWPC")  { input_path+="komen_singleM.dat"; }
  if (model=="SingleB")     { input_path+="komen_singleB.dat"; }
  if (model=="SingleC")     { input_path+="komen_singleC.dat"; }
  if (model=="SingleA")     { input_path+="komen_singleA.dat"; }
  if (model=="DoubleM1B_M") { input_path+="komen_doubleM1B_M.dat"; }
  if (model=="DoubleM1B_B") { input_path+="komen_doubleM1B_B.dat"; }
  if (model=="DoubleM2B_M") { input_path+="komen_doubleM2B_M.dat"; }
  if (model=="DoubleM2B_B") { input_path+="komen_doubleM2B_B.dat"; }
  if (model=="DoubleM3B_M") { input_path+="komen_doubleM3B_M.dat"; }
  if (model=="DoubleM3B_B") { input_path+="komen_doubleM3B_B.dat"; }
  if (model=="DoubleM4B_M") { input_path+="komen_doubleM4B_M.dat"; }
  if (model=="DoubleM4B_B") { input_path+="komen_doubleM4B_B.dat"; }
  if (model=="DoubleM5B_M") { input_path+="komen_doubleM5B_M.dat"; }
  if (model=="DoubleM5B_B") { input_path+="komen_doubleM5B_B.dat"; }
  if (model=="DoubleM1C_M") { input_path+="komen_doubleM1C_M.dat"; }
  if (model=="DoubleM1C_C") { input_path+="komen_doubleM1C_C.dat"; }
  if (model=="DoubleM2C_M") { input_path+="komen_doubleM2C_M.dat"; }
  if (model=="DoubleM2C_C") { input_path+="komen_doubleM2C_C.dat"; }
  if (model=="DoubleM3C_M") { input_path+="komen_doubleM3C_M.dat"; }
  if (model=="DoubleM3C_C") { input_path+="komen_doubleM3C_C.dat"; }
  if (model=="DoubleM4C_M") { input_path+="komen_doubleM4C_M.dat"; }
  if (model=="DoubleM4C_C") { input_path+="komen_doubleM4C_C.dat"; }
  if (model=="DoubleM5C_M") { input_path+="komen_doubleM5C_M.dat"; }
  if (model=="DoubleM5C_C") { input_path+="komen_doubleM5C_C.dat"; }
  if (model=="DoubleM1A_M") { input_path+="komen_doubleM1A_M.dat"; }
  if (model=="DoubleM1A_A") { input_path+="komen_doubleM1A_A.dat"; }
  if (model=="DoubleM2A_M") { input_path+="komen_doubleM2A_M.dat"; }
  if (model=="DoubleM2A_A") { input_path+="komen_doubleM2A_A.dat"; }
  if (model=="DoubleM3A_M") { input_path+="komen_doubleM3A_M.dat"; }
  if (model=="DoubleM3A_A") { input_path+="komen_doubleM3A_A.dat"; }
  if (model=="DoubleM4A_M") { input_path+="komen_doubleM4A_M.dat"; }
  if (model=="DoubleM4A_A") { input_path+="komen_doubleM4A_A.dat"; }
  if (model=="DoubleM5A_M") { input_path+="komen_doubleM5A_M.dat"; }
  if (model=="DoubleM5A_A") { input_path+="komen_doubleM5A_A.dat"; }
  if (model=="DoubleBC_B")  { input_path+="komen_doubleBC_B.dat"; }
  if (model=="DoubleBC_C")  { input_path+="komen_doubleBC_C.dat"; }
  if (model=="DoubleCA_C")  { input_path+="komen_doubleCA_C.dat"; }
  if (model=="DoubleCA_A")  { input_path+="komen_doubleCA_A.dat"; }
  if (model=="DoubleBC_BC") { input_path+="komen_doubleBC_B.dat"; }

  std::ifstream input_dat(input_path.c_str());
  if (!input_dat.is_open()) {
    SysError("Reconstruction", "Failed to open input file: %s", input_path.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  std::string input_path2 = output_dir + "/";
  std::ifstream input_dat2;
  if (model=="DoubleBC_BC") {
    input_path2+="komen_doubleBC_C.dat";
    input_dat2 = std::ifstream(input_path.c_str());
  }

  // 5. Initialize Reconstruction Tools
  SciBar100x100* detector = new SciBar100x100(0.0,0.0,0.0);
  detector -> SetLineColor(3);
  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  ClusterMaker<SciBar100x100*> cluster_maker(detector);
  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::CombinedLocal);
//  cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::CombinedGlobal);

  if (model=="DoubleBC_BC") {
    cluster_maker.SetClusterModel(ClusterMaker<SciBar100x100*>::Global);
  }

  SciBar100x100* detector2 = new SciBar100x100(0.0,0.0,0.0);
  TVector3 detector2_position;
  ClusterMaker<SciBar100x100*> cluster_maker2;
  if (model=="DoubleBC_BC") {
    detector2 -> MoveXYZ(detB_position.X()-1000.0,detB_position.Y(),detB_position.Z()+200.0);
    detector2 -> SetLineColor(3);
    detector2 -> SetFillColor(3);
    detector2 -> SetDensity(0.0);
    cluster_maker2 = ClusterMaker<SciBar100x100*>(detector2);
    cluster_maker2.SetClusterModel(ClusterMaker<SciBar100x100*>::Global);
  }


  TrackFinder track_finder;
  track_finder.SetTotalLayer(4);
  track_finder.SetChisqCut(1000000.0);
  track_finder.SetMinimumHits(4);
  if (model=="DoubleBC_BC") {
    track_finder.SetTotalLayer(16);
    track_finder.SetMinimumHits(8);
  }

//  TrackFinder track_finder(16, 1000000.0);
//  track_finder.SetMinimumHits(8);

//  TrackFinder track_finder(8, 1000000.0);
//  track_finder.SetMinimumHits(6);

//  TrackFinder track_finder(4, 1000000.0);
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

    std::vector<SpacePoint> space_points;
    for (const auto& space_point : all_space_points) {
      SpacePoint sp(space_point.X(), space_point.Y(), space_point.Z());
      sp.SetErrorX(space_point.GetErrorX());
      sp.SetErrorY(space_point.GetErrorY());
      sp.SetErrorZ(space_point.GetErrorZ());
      sp.SetUnitIndex(space_point.GetUnitIndex());
      space_points.push_back(sp);
    }

    if (model=="DoubleBC_BC") {
      std::getline(input_dat2, line);
      std::stringstream data2(line);
      std::vector<int> hit_ids2;
      while (data2 >> hit_id) { hit_ids2.push_back(hit_id); } 

      std::vector<SpacePoint> all_space_points2 = cluster_maker2.Execute(hit_ids2);
      for (const auto& space_point : all_space_points2) {
        SpacePoint sp(space_point.X(), space_point.Y(), space_point.Z());
        sp.SetErrorX(space_point.GetErrorX());
        sp.SetErrorY(space_point.GetErrorY());
        sp.SetErrorZ(space_point.GetErrorZ());
        sp.SetUnitIndex(space_point.GetUnitIndex()+8);
        space_points.push_back(sp);
      }
    }

    // Step B: Track Finding (Combinatorial Search and Fitting)
    Track best_track = track_finder.GetBestTrack(space_points);

// STSTST    Track best_track = track_finder.GetBestTrack(all_space_points);

    // Step C: Recording (Fill TTree if a valid track was found)
    if (best_track.IsValid()) {

      m_counter       = i_events;
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

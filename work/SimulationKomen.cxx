#ifdef __CLING__
R__ADD_INCLUDE_PATH(.);
#endif

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <TCanvas.h>
#include <TView.h>
#include <TColor.h>
#include <TError.h>
#include <TSystem.h>

#include "Simulation/ReadExpacs.h"
#include "Simulation/Muon.h"
#include "Simulation/MaterialList.h"
#include "Simulation/RectangularPrism.h"
#include "Simulation/RectangularPyramid.h"
#include "Simulation/Pipe.h"
#include "Simulation/House.h"
#include "Simulation/SciBar100x100.h"
#include "Geometry/Mountain.h"

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

//git add . → git commit -m "コメント" → git push 
//git pwd: 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
// source /home/takumu/root/bin/thisroot.sh
//root [0]  gSystem->AddIncludePath("-I${PWD}/../")
//root [1] .L SimulationKomen.cxx+
//root [2] SimulationKomen()

////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 PlaceDetector(Mountain* mountain,
                       double lat, double lon,
                       double ref_lat, double ref_lon,
                       double z_offset) {
    double x, y;
    mountain->reverse_explore(lat, lon, &x, &y, ref_lat, ref_lon);

    double z = mountain->Interpolate(x, y);

    return TVector3(x, y, z + z_offset);
}

int SimulationKomen(std::string model="ALL") {
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
  std::string pdffile      = "komen.pdf";

  double threshold = 1000.0;
  double geometry_rangexy = 50000.0;
  double geometry_rangez  = 2000.0;

  // Total Flux (no cut): mu+ : 0.008235 /cm2/s,   mu- : 0.007109 /cm2/s
  // Total Flux (E>1GeV): mu+ : 0.006446 /cm2/s,   mu- : 0.005470 /cm2/s
  double flux = 0.0;
  if (threshold==   0.0) { flux = 0.008235+0.007109; }
  if (threshold==1000.0) { flux = 0.006446+0.005470; }

  int n_events = flux*100.0*100.0*60*60*24;  // 1-day events
  if (model=="DoubleMB") {
    n_events = flux*100.0*100.0*60*60*1000/861576*24;  // 1-day events
  }
  if (model=="DoubleMA") {
    n_events = flux*100.0*100.0*60*60*1000/4889616*24;  // 1-day events
  }
  int n_draw = 20;

  double weight = 0.01;
  n_events *= weight;

  // TVector3 detA_position(0.0,0.0,0.0);
  // TVector3 detB_position(2000.0,0.0,0.0);
  // TVector3 MWPC_position(2000.0,0.0,2000.0);

  // 2. Setup Output
  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    ::Info("Simulation", "Created output directory: %s", output_dir.c_str());
  }

  // 3. Setup Geometry
  // RectangularPrism* ground = new RectangularPrism(1000.0, 0.0, 0.5*geometry_rangez, 
  //                                                 geometry_rangexy, geometry_rangexy, geometry_rangez);
  // ground -> SetLineColor(791);
  // ground -> SetFillColor(791);
  // ground -> SetDensity(GetDensity("NormalSoil"));
  double ref_lat = 31.9713;
  double ref_lon = 130.3420;

  Mountain* mountain = new Mountain();
  mountain->SetReference(ref_lat, ref_lon, 0.0);
  mountain->SetGeometry("../Geometry/Hokusatsu_elevation.dat");
  std::cout << "Geometry loaded" << std::endl;

  // A（地中）
  TVector3 detA_position = PlaceDetector(
      mountain, ref_lat, ref_lon, ref_lat, ref_lon, 250.0);

  // B（少しズラす）
  TVector3 detB_position = PlaceDetector(
      mountain, ref_lat, ref_lon, ref_lat, ref_lon, -50.0);

  // MWPC（上空）
  TVector3 MWPC_position = PlaceDetector(
      mountain, ref_lat, ref_lon, ref_lat, ref_lon, +1000.0);

      double x, y;
      mountain->reverse_explore(31.954, 130.324, &x, &y, ref_lat, ref_lon);

      std::cout << "x=" << x << " y=" << y << std::endl;

  

  

  SciBar100x100* detA = new SciBar100x100(0.0,0.0,0.0);
  detA -> MoveXYZ(detA_position.X(),detA_position.Y(),detA_position.Z());
  detA -> SetLineColor(3);
  detA -> SetFillColor(3);
  detA -> SetDensity(0.0);

  SciBar100x100* detB = new SciBar100x100(0.0,0.0,0.0);
  detB -> MoveXYZ(detB_position.X(),detB_position.Y(),detB_position.Z());
//  detB -> RotateTheta(TMath::Pi()/3.0);
//  detB -> RotatePhi(TMath::Pi()/4.0);
  detB -> SetLineColor(3);
  detB -> SetFillColor(3);
  detB -> SetDensity(0.0);

  SciBar100x100* MWPC = new SciBar100x100(0.0,0.0,0.0);
  MWPC -> MoveXYZ(MWPC_position.X(),MWPC_position.Y(),MWPC_position.Z());
  MWPC -> SetLineColor(3);
  MWPC -> SetFillColor(3);
  MWPC -> SetDensity(0.0);

  MWPC   -> SetPriority(1);
  detB   -> SetPriority(2);
  detA   -> SetPriority(3);
  // ground -> SetPriority(4);
  std::vector<Shape*> shapes;

 // shapes.push_back(ground);
  shapes.push_back(MWPC);
  shapes.push_back(detB);
  shapes.push_back(detA);

  // 4. Setup Muon Generator
  ReadExpacs* reader = new ReadExpacs(threshold);

  int n = std::to_string(n_events).length();
  int n_generation_failed = 0;
  int n_below_threshold = 0;
  int i_draw = 0;

  int n_generation_success = 0;
  int n_total_attempts     = 0;
  int n_above_threshold    = 0;
  int n_nohit_detector     = 0;
  int n_hit_detector       = 0;

  // Monitoring tracks
  std::vector<Muon*> muonA_tracks;
  std::vector<Muon*> muonB_tracks;
  std::vector<Muon*> muonM_tracks;
  std::vector<Muon*> muonMB_tracks;
  std::vector<Muon*> muonMA_tracks;

  //==============================================================================
  // Simulate detector A
  //==============================================================================
  std::string outputA_path = output_dir + "/" + datafileA; 
  std::ofstream outputA_dat(outputA_path); 

  if (model=="ALL" || model=="SingleA") { 
    if (!outputA_dat.is_open()) { 
        SysError("Simulation::Initialize", "Could not open output file: %s", datafileA.c_str()); 
        gSystem->Exit(EXIT_FAILURE); 
    }

    for (int i = 0; i < n_events; ++i) {
        Muon* mu = new Muon(reader, 12345);
        int trial = 0;
        while (!mu->Generate() && trial < 10000) {
          trial++;
        }
        auto hits = detA->IsHitList(*mu);
        if (!hits.empty()) {
          for (int id : hits) {
          outputA_dat << id << " ";
       }
       outputA_dat << "\n";
      }
      
    }
    outputA_dat.close();
  }
    
    Muon* muonA = new Muon(reader, 123456);
    muonA->SetReferenceZ(detA_position.Z()+15.0);
    muonA->SetEndZ(-500.0);
    muonA->SetStartZ(3000.0);
    muonA->Set1stBoundaryXY(detA_position.X()-50.0, detA_position.X()+50.0,
                            detA_position.Y()-50.0, detA_position.Y()+50.0);

    // 5. Event Loop
    n = std::to_string(n_events).length();
    n_generation_failed = 0;
    n_below_threshold = 0;
    i_draw = 0;
    for (int i = 0; i < n_events; ++i) {
      if (i%(n_events/1000)==0) {         
        Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events); 
      }

      int trial = 0;
      while (!muonA->Generate() && trial < 10000) {
          n_generation_failed++;
          trial++;
      }

      auto [initial_energy, loss_energy] = muonA->CalculateEnergyLoss(shapes);
      double final_energy = initial_energy - loss_energy;
      if (g_debug_level>0) {
        Info("Simulation", "Muon energy (initial, loss, final): (%.2f, %.2f, %.2f)", 
            initial_energy, loss_energy, final_energy);
      }

      if (final_energy <= threshold) {
        n_below_threshold++;
        continue;
      }

      // Store hits in detector
      std::vector<int> all_hit_ids = detA->IsHitList(*muonA);
      for (int id : all_hit_ids) {
        outputA_dat << std::setw(5) << id;
      }
      outputA_dat << "\n";      

      if (i_draw < n_draw) {
        Muon* track = new Muon(*muonA);
        track -> SetMarkerColor(kPink);
        track -> SetLineColor(kViolet);
        detA -> MarkIntersections(track);
        muonA_tracks.push_back(track);
        i_draw++;
      }
    }
    Info("Simulation::EventLoop", "Progress: %*d / %d completed", n, n_events, n_events);

    // 6. Final Statistics
    n_generation_success = n_events;
    n_total_attempts     = n_generation_success + n_generation_failed;
    n_above_threshold    = n_generation_success - n_below_threshold;

    Info("Simulation::Statistics", "[Generation] Attempts : %*d", n, n_total_attempts);
    Info("Simulation::Statistics", "[Generation] Discarded: %*d", n, n_generation_failed);
    Info("Simulation::Statistics", "[Generation] Generated: %*d", n, n_generation_success);
    Info("Simulation::Statistics", "[Threshold]  Filtered : %*d", n, n_below_threshold);
    Info("Simulation::Statistics", "[Threshold]  Passed   : %*d", n, n_above_threshold);

    // 7. Finalize
    outputA_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputA_path.c_str());

  //==============================================================================
  // Simulate detector B
  //==============================================================================
  std::string outputB_path = output_dir + "/" + datafileB; 
  std::ofstream outputB_dat(outputB_path); 
  if (model=="ALL" || model=="SingleB") { 
    if (!outputB_dat.is_open()) { 
        SysError("Simulation::Initialize", "Could not open output file: %s", datafileB.c_str()); 
        gSystem->Exit(EXIT_FAILURE); 
    }
    for (int i = 0; i < n_events; ++i) {
        Muon* mu = new Muon(reader,12345);
        int trial = 0;
        while (!mu->Generate() && trial < 10000) {
            trial++;
        }
        auto hits = detA->IsHitList(*mu);
        if (!hits.empty()) {
          for (int id : hits) {
              outputB_dat << id << " ";
          }
        outputB_dat << "\n";
        }
    }
    outputB_dat.close();
  }

    Muon* muonB = new Muon(reader, 7645311);
    muonB->SetReferenceZ(detB_position.Z()+15.0);
    muonB->SetEndZ(-500.0);
    muonB->SetStartZ(3000.0);
    muonB->Set1stBoundaryXY(detB_position.X()-50.0, detB_position.X()+50.0,
                            detB_position.Y()-50.0, detB_position.Y()+50.0);

    // 5. Event Loop
    n_generation_failed = 0;
    n_below_threshold = 0;
    i_draw = 0;
    for (int i = 0; i < n_events; ++i) {
      if (i%(n_events/10)==0) {
        Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);
      }

      int trial = 0;
      while (!muonB->Generate() && trial < 10000) {
          n_generation_failed++;
          trial++;
      }

      auto [initial_energy, loss_energy] = muonB->CalculateEnergyLoss(shapes);
      double final_energy = initial_energy - loss_energy;
      if (g_debug_level>0) {
        Info("Simulation", "Muon energy (initial, loss, final): (%.2f, %.2f, %.2f)", 
            initial_energy, loss_energy, final_energy);
      }

      if (final_energy <= threshold) {
        n_below_threshold++;
        continue;
      }

      // Store hits in detector
      std::vector<int> all_hit_ids = detB->IsHitList(*muonB);
      for (int id : all_hit_ids) {
        outputB_dat << std::setw(5) << id;
      }
      outputB_dat << "\n";      

      if (i_draw < n_draw) {
        Muon* track = new Muon(*muonB);
        track -> SetMarkerColor(kPink);
        track -> SetLineColor(kViolet);
        detB -> MarkIntersections(track);
        muonB_tracks.push_back(track);
        i_draw++;
      }
    }
    Info("Simulation::EventLoop", "Progress: %*d / %d completed", n, n_events, n_events);

    // 6. Final Statistics
    n_generation_success = n_events;
    n_total_attempts     = n_generation_success + n_generation_failed;
    n_above_threshold    = n_generation_success - n_below_threshold;

    Info("Simulation::Statistics", "[Generation] Attempts : %*d", n, n_total_attempts);
    Info("Simulation::Statistics", "[Generation] Discarded: %*d", n, n_generation_failed);
    Info("Simulation::Statistics", "[Generation] Generated: %*d", n, n_generation_success);
    Info("Simulation::Statistics", "[Threshold]  Filtered : %*d", n, n_below_threshold);
    Info("Simulation::Statistics", "[Threshold]  Passed   : %*d", n, n_above_threshold);

    // 7. Finalize
    outputB_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputB_path.c_str());

  //==============================================================================
  // Simulate detector MWPC
  //==============================================================================
  std::string outputM_path = output_dir + "/" + datafileM; 
  std::ofstream outputM_dat(outputM_path); 
  if (model=="ALL" || model=="SingleMWPC") { 
    if (!outputM_dat.is_open()) { 
        SysError("Simulation::Initialize", "Could not open output file: %s", datafileM.c_str()); 
        gSystem->Exit(EXIT_FAILURE); 
    }
    for (int i = 0; i < n_events; ++i) {
        Muon* mu = new Muon(reader, 12345);
        int trial = 0;
        while (!mu->Generate() && trial < 10000) {
            trial++;
        }
        auto hits = detA->IsHitList(*mu);
        if (!hits.empty()) {
          for (int id : hits) {
            outputM_dat << id << " ";
          }
          outputM_dat << "\n";
        }
    }
    outputM_dat.close();
  }

    Muon* muonM = new Muon(reader, 54738291);
    muonM->SetReferenceZ(MWPC_position.Z()+15.0);
    muonM->SetEndZ(-500.0);
    muonM->SetStartZ(3000.0);
    muonM->Set1stBoundaryXY(MWPC_position.X()-50.0, MWPC_position.X()+50.0,
                            MWPC_position.Y()-50.0, MWPC_position.Y()+50.0);

    // 5. Event Loop
    n_generation_failed = 0;
    n_below_threshold = 0;
    i_draw = 0;
    for (int i = 0; i < n_events; ++i) {
      if (i%(n_events/10)==0) {
        Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);
      }

    int trial = 0;
    while (!muonM->Generate() && trial < 10000) {
        n_generation_failed++;
        trial++;
    }

      // Store hits in detector
      std::vector<int> all_hit_ids = MWPC->IsHitList(*muonM);
      for (int id : all_hit_ids) {
        outputM_dat << std::setw(5) << id;
      }
      outputM_dat << "\n";      

      if (i_draw < n_draw) {
        Muon* track = new Muon(*muonM);
        track -> SetMarkerColor(kPink);
        track -> SetLineColor(kViolet);
        MWPC -> MarkIntersections(track);
        muonM_tracks.push_back(track);
        i_draw++;
      }
    }
    Info("Simulation::EventLoop", "Progress: %*d / %d completed", n, n_events, n_events);

    // 6. Final Statistics
    n_generation_success = n_events;
    n_total_attempts     = n_generation_success + n_generation_failed;
    n_above_threshold    = n_generation_success - n_below_threshold;

    Info("Simulation::Statistics", "[Generation] Attempts : %*d", n, n_total_attempts);
    Info("Simulation::Statistics", "[Generation] Discarded: %*d", n, n_generation_failed);
    Info("Simulation::Statistics", "[Generation] Generated: %*d", n, n_generation_success);
    Info("Simulation::Statistics", "[Threshold]  Filtered : %*d", n, n_below_threshold);
    Info("Simulation::Statistics", "[Threshold]  Passed   : %*d", n, n_above_threshold);

    // 7. Finalize
    outputM_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputM_path.c_str());
  

  //==============================================================================
  // Simulate simultaneous hits in MWPC and detector B
  //==============================================================================
  std::string outputMB_M_path = output_dir + "/" + datafileMB_M; 
  std::ofstream outputMB_M_dat(outputMB_M_path); 
  if (model=="ALL" || model=="DoubleMB") { 
    if (!outputMB_M_dat.is_open()) { 
        SysError("Simulation::Initialize", "Could not open output file: %s", datafileMB_M.c_str()); 
        gSystem->Exit(EXIT_FAILURE); 
    }
    for (int i = 0; i < n_events; ++i) {
        Muon* mu = new Muon(reader, 12345);
        int trial = 0;
        while (!mu->Generate() && trial < 10000) {
            trial++;
        }
        auto hits = MWPC->IsHitList(*mu);
        if (!hits.empty()) {
          for (int id : hits) {
          outputMB_M_dat << id << " ";
        }
        outputMB_M_dat << "\n";
      }
    }
    outputMB_M_dat.close();
  }

    std::string outputMB_B_path = output_dir + "/" + datafileMB_B;
    std::ofstream outputMB_B_dat(outputMB_B_path);
    if (!outputMB_B_dat.is_open()) {
      SysError("Simulation::Initialize", "Could not open output file: %s", datafileMB_B.c_str());
      gSystem->Exit(EXIT_FAILURE);
    }

    Muon* muonMB = new Muon(reader, 9918291);
    muonMB->SetReferenceZ(detB_position.Z()+15.0);
    muonMB->SetEndZ(-500.0);
    muonMB->SetStartZ(3000.0);
    muonMB->Set1stBoundaryXY(detB_position.X()-50.0, detB_position.X()+50.0,
        detB_position.Y()-50.0, detB_position.Y()+50.0);

    muonMB->SetSecondReferenceZ(MWPC_position.Z()+15.0);
    muonMB->Set2ndBoundaryXY(MWPC_position.X()-50.0, MWPC_position.X()+50.0,
                             MWPC_position.Y()-50.0, MWPC_position.Y()+50.0);

    // 5. Event Loop
    n_generation_failed = 0;
    n_below_threshold = 0;
    n_nohit_detector = 0;
    i_draw = 0;
    for (int i = 0; i < n_events; ++i) {
      if (i % (n_events / 10) == 0) Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);

      int trial = 0;
      while (!muonMB->Generate() && trial < 10000) {
          n_generation_failed++;
          trial++;
      }

      auto [initial_energy, loss_energy] = muonMB->CalculateEnergyLoss(shapes);
      double final_energy = initial_energy - loss_energy;
      if (g_debug_level>0) {
        Info("Simulation", "Muon energy (initial, loss, final): (%.2f, %.2f, %.2f)", initial_energy, loss_energy, final_energy);
      }

      if (final_energy <= threshold) {
        n_below_threshold++;
        continue;
      }

      // Store hits in detector
      std::vector<int> mwpc_hit_ids = MWPC->IsHitList(*muonMB);
      std::vector<int> detB_hit_ids = detB->IsHitList(*muonMB);
      if (mwpc_hit_ids.size()==0 || detB_hit_ids.size()==0) {
        n_nohit_detector++;
        continue;
      }

      for (int id : mwpc_hit_ids) {
        outputMB_M_dat << std::setw(5) << id;
      }
      outputMB_M_dat << "\n";      

      for (int id : detB_hit_ids) {
        outputMB_B_dat << std::setw(5) << id;
      }
      outputMB_B_dat << "\n";      

      if (i_draw < n_draw) {
        Muon* track = new Muon(*muonMB);
        track -> SetMarkerColor(kPink);
        track -> SetLineColor(4);
        muonMB_tracks.push_back(track);
        i_draw++;
      }
    }
    Info("Simulation::EventLoop", "Progress: %*d / %d completed", n, n_events, n_events);

    // 6. Final Statistics
    n_generation_success = n_events;
    n_total_attempts     = n_generation_success + n_generation_failed;
    n_above_threshold    = n_generation_success - n_below_threshold;
    n_hit_detector       = n_above_threshold - n_nohit_detector;

    Info("Simulation::Statistics", "[Generation]  Attempts : %*d", n, n_total_attempts);
    Info("Simulation::Statistics", "[Generation]  Discarded: %*d", n, n_generation_failed);
    Info("Simulation::Statistics", "[Generation]  Generated: %*d", n, n_generation_success);
    Info("Simulation::Statistics", "[Threshold]   Filtered : %*d", n, n_below_threshold);
    Info("Simulation::Statistics", "[Threshold]   Passed   : %*d", n, n_above_threshold);
    Info("Simulation::Statistics", "[DetectorHit] Filtered : %*d", n, n_nohit_detector);
    Info("Simulation::Statistics", "[DetectorHit] Passed   : %*d", n, n_hit_detector);

    // 7. Finalize
    outputMB_M_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputMB_M_path.c_str());

    outputMB_B_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputMB_B_path.c_str());
  

  //==============================================================================
  // Simulate simultaneous hits in MWPC and detector A
  //==============================================================================
  std::string outputMA_M_path = output_dir + "/" + datafileMA_M; 
  std::ofstream outputMA_M_dat(outputMA_M_path);
  if (model=="ALL" || model=="DoubleMA") {  
    if (!outputMA_M_dat.is_open()) { 
        SysError("Simulation::Initialize", "Could not open output file: %s", datafileMA_M.c_str()); 
        gSystem->Exit(EXIT_FAILURE); 
    }
    for (int i = 0; i < n_events; ++i) {
        Muon* mu = new Muon(reader, 12345);
        int trial = 0;
        while (!mu->Generate() && trial < 10000) {
            trial++;
        }
        auto hitsM = MWPC->IsHitList(*mu);
        auto hitsA = detA->IsHitList(*mu);

        if (!hitsM.empty() && !hitsA.empty()) {
          for (int id : hitsM) {
          outputMA_M_dat << id << " ";
        }
      outputMA_M_dat << "\n";
      }
    }
    outputMA_M_dat.close();
  }

    std::string outputMA_A_path = output_dir + "/" + datafileMA_A;
    std::ofstream outputMA_A_dat(outputMA_A_path);
    if (!outputMA_A_dat.is_open()) {
      SysError("Simulation::Initialize", "Could not open output file: %s", datafileMA_A.c_str());
      gSystem->Exit(EXIT_FAILURE);
    }

    Muon* muonMA = new Muon(reader, 39781);
    muonMA->SetReferenceZ(detA_position.Z()+15.0);
    muonMA->SetEndZ(-500.0);
    muonMA->SetStartZ(3000.0);
    muonMA->Set1stBoundaryXY(detA_position.X()-50.0, detA_position.X()+50.0,
        detA_position.Y()-50.0, detA_position.Y()+50.0);

    muonMA->SetSecondReferenceZ(MWPC_position.Z()+15.0);
    muonMA->Set2ndBoundaryXY(MWPC_position.X()-50.0, MWPC_position.X()+50.0,
        MWPC_position.Y()-50.0, MWPC_position.Y()+50.0);


    // 5. Event Loop
    n_generation_failed = 0;
    n_below_threshold = 0;
    n_nohit_detector = 0;
    i_draw = 0;
    for (int i = 0; i < n_events; ++i) {
      if (i%(n_events/10)==0) {
        Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);
      }

      int trial = 0;
      while (!muonMA->Generate() && trial < 10000) {
          n_generation_failed++;
          trial++;
      }

      auto [initial_energy, loss_energy] = muonMA->CalculateEnergyLoss(shapes);
      double final_energy = initial_energy - loss_energy;
      if (g_debug_level>0) {
        Info("Simulation", "Muon energy (initial, loss, final): (%.2f, %.2f, %.2f)", initial_energy, loss_energy, final_energy);
      }

      if (final_energy <= threshold) {
        n_below_threshold++;
        continue;
      }

      // Store hits in detector
      std::vector<int> mwpc_hit_ids = MWPC->IsHitList(*muonMA);
      std::vector<int> detA_hit_ids = detA->IsHitList(*muonMA);
      if (mwpc_hit_ids.size()==0 || detA_hit_ids.size()==0) {
        n_nohit_detector++;
        continue;
      }

      for (int id : mwpc_hit_ids) {
        outputMA_M_dat << std::setw(5) << id;
      }
      outputMA_M_dat << "\n";      

      for (int id : detA_hit_ids) {
        outputMA_A_dat << std::setw(5) << id;
      }
      outputMA_A_dat << "\n";      

      if (i_draw < n_draw) {
        Muon* track = new Muon(*muonMA);
        track -> SetMarkerColor(kPink);
        track -> SetLineColor(4);
        muonMA_tracks.push_back(track);
        i_draw++;
      }
    }
    Info("Simulation::EventLoop", "Progress: %*d / %d completed", n, n_events, n_events);

    // 6. Final Statistics
    n_generation_success = n_events;
    n_total_attempts     = n_generation_success + n_generation_failed;
    n_above_threshold    = n_generation_success - n_below_threshold;
    n_hit_detector       = n_above_threshold - n_nohit_detector;

    Info("Simulation::Statistics", "[Generation]  Attempts : %*d", n, n_total_attempts);
    Info("Simulation::Statistics", "[Generation]  Discarded: %*d", n, n_generation_failed);
    Info("Simulation::Statistics", "[Generation]  Generated: %*d", n, n_generation_success);
    Info("Simulation::Statistics", "[Threshold]   Filtered : %*d", n, n_below_threshold);
    Info("Simulation::Statistics", "[Threshold]   Passed   : %*d", n, n_above_threshold);
    Info("Simulation::Statistics", "[DetectorHit] Filtered : %*d", n, n_nohit_detector);
    Info("Simulation::Statistics", "[DetectorHit] Passed   : %*d", n, n_hit_detector);

    // 7. Finalize
    outputMA_M_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputMA_M_path.c_str());

    outputMA_A_dat.close();
    Info("Simulation::Finalize", "Output File: %s", outputMA_A_path.c_str());

  //==============================================================================

  // 8. Visualize
  int irep;
  TCanvas* canvas = new TCanvas("canvas", "Simulation", 0, 0, 800, 800);
  TView* view = TView::CreateView();
  view->SetRange(-20000.0, -20000.0, 0.0, 20000.0, 20000.0, 30000.0);
  view->SetView(800.0, 800.0, 0.0, irep);
  mountain -> DrawWithMap();
  

  for (auto shape : shapes) shape->Draw("f");
  if (model=="ALL" || model=="SingleA") {
    for (auto track : muonA_tracks) { track->Draw("m"); }
  }
  if (model=="ALL" || model=="SingleB") {
    for (auto track : muonB_tracks) { track->Draw("m"); }
  }
  if (model=="ALL" || model=="SingleMWPC") {
    for (auto track : muonM_tracks) { track->Draw("m"); }
  }
  if (model=="ALL" || model=="DoubleMB") {
    for (auto track : muonMB_tracks) { track->Draw("m"); }
  }
  if (model=="ALL" || model=="DoubleMA") {
    for (auto track : muonMA_tracks) { track->Draw("m"); }
  }

  canvas->Modified();
  canvas->Update();
  std::string output_path = output_dir + "/" + pdffile;
  canvas->SaveAs(output_path.c_str());

  return 0;
}

#ifndef __CLING__
int main(int argc, char** argv) {
  std::string model = argv[1];
  return SimulationKomen(model);
}
#endif

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
#include "Simulation/TriangularPrism.h"
#include "Simulation/Pipe.h"
#include "Simulation/SciBar100x100.h"

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

int SimulationTest(bool withAnomaly=true) {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

  std::string output_dir = "data/output/test";

//  bool withAnomaly = true;

  std::string datafile = "test_anomaly.dat";
  std::string pdffile  = "test_anomaly.pdf";
  int randomseed = 1234567;
  if (!withAnomaly) {
    datafile = "test_noanomaly.dat";
    pdffile  = "test_noanomaly.pdf";
    randomseed = 762149;
  }

//  int n_events = 1000000;
  int n_events = 500000;
  int n_draw = 20;

  double threshold = 110.0;

  // 2. Setup Output
  std::string output_path = output_dir + "/" + datafile;
  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    ::Info("Simulation", "Created output directory: %s", output_dir.c_str());
  }
  std::ofstream output_dat(output_path);
  if (!output_dat.is_open()) {
    SysError("Simulation::Initialize", "Could not open output file: %s", datafile.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  // 3. Setup Geometry
//  TriangularPrism* anomaly = new TriangularPrism(  0.0,   0.0, 50.0, 
//                                                  40.0,  20.0*TMath::Sqrt(3.0), 50.0, false);
  TriangularPrism* anomaly = new TriangularPrism(  0.0,   0.0, 75.0, 
                                                  40.0,  20.0*TMath::Sqrt(3.0), 100.0, false);
//  anomaly -> RotateTheta(TMath::Pi());
//  anomaly -> MoveXYZ(100.0, 0.0, 0.0);
  anomaly -> SetLineColor(kAzure-4);
  anomaly -> SetFillColor(kAzure-4);
  anomaly -> SetMaterial(MaterialList::GetMaterial("Steel"));
//  anomaly -> SetDensity(0.0);

  SciBar100x100* detector = new SciBar100x100(0.0,0.0,0.0);
  detector -> SetLineColor(3);
  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  anomaly  -> SetPriority(2);
  detector -> SetPriority(1);
  std::vector<Shape*> shapes;

  shapes.push_back(detector);
  if (withAnomaly) {
    shapes.push_back(anomaly);
  }

  // 4. Setup Muon Generator
  ReadExpacs* reader = new ReadExpacs(0.0);

  Muon* muon = new Muon(reader, randomseed);
  muon->SetReferenceZ(15.0);
  muon->SetEndZ(-100.0);
  muon->SetStartZ(500.0);
  muon->Set1stBoundaryXY(-50.0, 50.0, -50.0, 50.0);

  // 5. Event Loop
  std::vector<Muon*> muon_tracks;
  int n = std::to_string(n_events).length();
  int n_generation_failed = 0;
  int n_below_threshold = 0;
  int i_draw = 0;
  for (int i = 0; i < n_events; ++i) {
    if (i % (n_events / 10) == 0) Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);

//    while (!muon->Generate()) {
//      n_generation_failed++;
//    }

    for (;;) {
      if (!muon->Generate() || muon->GetEnergyMuon()>2000.0) {
        n_generation_failed++;
      }
      else {
        break;
      }
    }

    auto [initial_energy, loss_energy] = muon->CalculateEnergyLoss(shapes);
    double final_energy = initial_energy - loss_energy;
    if (g_debug_level>0) {
      Info("Simulation", "Muon energy (initial, loss, final): (%.2f, %.2f, %.2f)", initial_energy, loss_energy, final_energy);
    }

    if (final_energy <= threshold) {
      n_below_threshold++;
      continue;
    }

//    if (loss_energy>0.0) {
//      std::cout << "Energy loss " << loss_energy << " " << initial_energy << " " << final_energy << std::endl;
//    }

    // Store hits in detector
    std::vector<int> all_hit_ids = detector->IsHitList(*muon);
    output_dat << std::fixed << std::setprecision(4) 
               << std::setw(12) << initial_energy 
               << std::setw(12) << loss_energy
               << std::setw(3) << muon->GetCharge() 
               << std::setw(8) << muon->GetTheta() 
               << std::setw(8) << muon->GetPhi();
    for (int id : all_hit_ids) {
      output_dat << std::setw(5) << id;
    }
    output_dat << "\n";      

    if (i_draw < n_draw) {
      Muon* track = new Muon(*muon);
      track -> SetMarkerColor(kPink);
      track -> SetLineColor(kViolet);
      detector -> MarkIntersections(track);
      muon_tracks.push_back(track);
      i_draw++;
    }
  }
  Info("Simulation::EventLoop", "Progress: %*d / %d completed", n, n_events, n_events);

  // 6. Final Statistics
  int n_generation_success = n_events;
  int n_total_attempts     = n_generation_success + n_generation_failed;
  int n_above_threshold    = n_generation_success - n_below_threshold;

  Info("Simulation::Statistics", "[Generation] Attempts : %*d", n, n_total_attempts);
  Info("Simulation::Statistics", "[Generation] Discarded: %*d", n, n_generation_failed);
  Info("Simulation::Statistics", "[Generation] Generated: %*d", n, n_generation_success);
  Info("Simulation::Statistics", "[Threshold]  Filtered : %*d", n, n_below_threshold);
  Info("Simulation::Statistics", "[Threshold]  Passed   : %*d", n, n_above_threshold);

  // 7. Finalize
  output_dat.close();
  Info("Simulation::Finalize", "Output File: %s", output_path.c_str());

  // 8. Visualize
  int irep;
  TCanvas* canvas = new TCanvas("canvas", "Simulation", 0, 0, 800, 800);
  TView* view = TView::CreateView();
  view->SetRange(-200.0, -200.0, 0.0, 200.0, 200.0, 400.0);
  view->SetView(200.0, 200.0, 0.0, irep);

  for (auto shape : shapes) shape->Draw("f");
  for (auto track : muon_tracks) track->Draw("m");

  canvas->Modified();
  canvas->Update();
  output_path = output_dir + "/" + pdffile;
  canvas->SaveAs(output_path.c_str());

  return 0;
}

#ifndef __CLING__
int main(int argc, char** argv) {
  std::string chk = argv[1];
  bool withAnomaly=true;
  if (chk.find("false")!=std::string::npos) { withAnomaly=false; }
  return SimulationTest(withAnomaly);
}
#endif

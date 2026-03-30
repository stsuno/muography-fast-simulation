#ifdef __CLING__
R__ADD_INCLUDE_PATH(.);
#endif

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

#include <TCanvas.h>
#include <TView.h>
#include <TColor.h>
#include <TError.h>
#include <TSystem.h>

#include "Simulation/ReadExpacs.h"
#include "Simulation/Muon.h"
#include "Simulation/MaterialList.h"

#include "Simulation/SciFiberCylinder.h"
#include "Simulation/SciBar50x50.h"

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

int SimulationFiberDetector() {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

  std::string output_dir = "data/output/fiberdetector";

  std::string datafile = "fiberdetector.dat";
  std::string pdffile  = "fiberdetector.pdf";

  int n_events = 1000;
  int n_draw = 100;

  double threshold = 0.0;
  double geometry_rangexy = 20.0;

  TVector3 detector_position(0.0,0.0,0.0);

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

  SciFiberCylinder* detector = new SciFiberCylinder(0.0, 0.0, 0.0);
//  detector -> MoveXYZ(detector_position.X(),detector_position.Y(),detector_position.Z());
//  detector -> RotateTheta(TMath::Pi()/3.0);
//  detector -> RotatePhi(TMath::Pi()/4.0);
//  detector -> SetLineColor(3);
//  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  detector -> SetPriority(1);
  std::vector<Shape*> shapes;

  shapes.push_back(detector);

  // 4. Setup Muon Generator
  ReadExpacs* reader = new ReadExpacs(threshold);

  Muon* muon = new Muon(reader);
//  muon->SetReferenceZ(detector_position.Z()-15.0);
  muon->SetReferenceZ(detector_position.Z());
//  muon->SetEndZ(-50.0);
//  muon->SetStartZ(200.0);
//  muon->Set1stBoundaryXY(detector_position.X()-25.0, detector_position.X()+25.0,
//                         detector_position.Y()-25.0, detector_position.Y()+25.0);

  double height = detector->GetDetectorHeight();
  muon -> SetCylinderBoundary(10.0, height, 1);    // outer boundary
//  muon->Set1stBoundaryXY(detector_position.X()-5.0, detector_position.X()+5.0,
//                         detector_position.Y()-5.0, detector_position.Y()+5.0);

  muon->SetStartZ(2.0*height);
  muon->SetEndZ(-2.0*height);

  // 5. Event Loop
  std::vector<Muon*> muon_tracks;
  int n = std::to_string(n_events).length();
  int n_generation_failed = 0;
  int n_below_threshold = 0;
  int n_hit_detector = 0;
  int i_draw = 0;
  for (int i = 0; i < n_events; ++i) {
    if (i % (n_events / 10) == 0) Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);

    while (!muon->Generate()) { n_generation_failed++; }

    auto [initial_energy, loss_energy] = muon->CalculateEnergyLoss(shapes);
    double final_energy = initial_energy - loss_energy;
    if (g_debug_level>0) {
      Info("Simulation", "Muon energy (initial, loss, final): (%.2f, %.2f, %.2f)", initial_energy, loss_energy, final_energy);
    }

    if (final_energy <= threshold) {
      n_below_threshold++;
      continue;
    }

    // Store hits in detector
    std::vector<std::pair<int, std::vector<int>>> all_hit_ids = detector->IsHitList(*muon);
    if (all_hit_ids.size()>0) {
      n_hit_detector++;
      double element = detector->GetFiber(0)->GetElementLength();
      double speed   = detector->GetFiber(0)->GetSignalSpeed();
      for (auto id : all_hit_ids) {
        int minID = detector->GetFiber(id.first)->GetHitID(*muon);

        TVector3 hitCenter(0.0,0.0,0.0);
        if (minID>0) { hitCenter += detector->GetFiber(id.first)->GetHitIDCenter(minID-1); }

        TVector3 hitPosition = detector->GetFiber(id.first)->GetHitPosition(*muon);
        double length  = element*(minID-1.5)+(hitCenter-hitPosition).Mag();
        output_dat << std::setw(5) << id.first << std::setw(10) << length/speed;
      }
      output_dat << "\n";      
    }

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
  Info("Simulation::Statistics", "[Detector]   Hits     : %*d", n, n_hit_detector);

  // 7. Finalize
  output_dat.close();
  Info("Simulation::Finalize", "Output File: %s", output_path.c_str());

  // 8. Visualize
  int irep;
  TCanvas* canvas = new TCanvas("canvas", "Simulation", 0, 0, 800, 800);
  TView* view = TView::CreateView();
  view->SetRange(-geometry_rangexy/2, -geometry_rangexy/2, -25.0, geometry_rangexy/2, geometry_rangexy/2, 25.0);
//  view->SetView(800.0, 800.0, 0.0, irep);
  view->SetView(50.0, 50.0, 0.0, irep);

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
  return SimulationFiberDetector();
}
#endif

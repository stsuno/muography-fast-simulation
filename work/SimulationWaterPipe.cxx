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
#include "Simulation/RectangularPrism.h"
#include "Simulation/RectangularPyramid.h"
#include "Simulation/Pipe.h"
#include "Simulation/House.h"
#include "Simulation/SciBar50x50.h"

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

int SimulationWaterPipe() {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

  std::string output_dir = "data/output/underground";

  std::string datafile = "underground.dat";
  std::string pdffile  = "underground.pdf";

  int n_events = 1000;
  int n_draw = 20;

  double threshold = 0.0;
  double geometry_rangexy = 1600.0;

  TVector3 detector_position(100.0,-220.0,15.0);

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
  RectangularPrism* ground = new RectangularPrism(0.0, 0.0, 150.0, geometry_rangexy, geometry_rangexy, 300.0);
  ground -> SetLineColor(791);
  ground -> SetFillColor(791);
  ground -> SetDensity(GetDensity("NormalSoil"));
  ground -> SetAtomicNumber(GetAtomicNumber("StandardRock"));
  ground -> SetAtomicMass(GetAtomicMass("StandardRock"));
  ground -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("StandardRock"));
  ground -> SetSternheimerParameters(GetSternheimerParameters("StandardRock"));

  RectangularPrism* road = new RectangularPrism(0.0, 0.0, 305.0, geometry_rangexy, geometry_rangexy, 10.0);
  road -> SetLineColor(kGray+1);
  road -> SetFillColor(kGray);
  road -> SetDensity(GetDensity("ShieldingConcrete"));
  road -> SetAtomicNumber(GetAtomicNumber("ShieldingConcrete"));
  road -> SetAtomicMass(GetAtomicMass("ShieldingConcrete"));
  road -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("ShieldingConcrete"));
  road -> SetSternheimerParameters(GetSternheimerParameters("ShieldingConcrete"));

  RectangularPrism* road2 = new RectangularPrism(0.0, 200.0, 310.0, geometry_rangexy, 20.0, 20.0);
  road2 -> SetLineColor(kGray+1);
  road2 -> SetFillColor(kGray);
  road2 -> SetDensity(GetDensity("ShieldingConcrete"));
  road2 -> SetAtomicNumber(GetAtomicNumber("ShieldingConcrete"));
  road2 -> SetAtomicMass(GetAtomicMass("ShieldingConcrete"));
  road2 -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("ShieldingConcrete"));
  road2 -> SetSternheimerParameters(GetSternheimerParameters("ShieldingConcrete"));

  House *house = new House(200.0, 600.0, 310.0, 800.0, 400.0, 300.0);
  house -> SetDensity(GetDensity("StandardRock"));
  house -> SetAtomicNumber(GetAtomicNumber("StandardRock"));
  house -> SetAtomicMass(GetAtomicMass("StandardRock"));
  house -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("StandardRock"));
  house -> SetSternheimerParameters(GetSternheimerParameters("StandardRock"));
  house -> SetRoofColor(901);

  Pipe* pipe = new Pipe(0.0, 0.0, 200.0, geometry_rangexy, 14.0, 15.0);
  pipe -> RotateTheta(TMath::Pi()/2.0);
  pipe -> RotatePhi(TMath::Pi()/2.0);
  pipe -> SetLineColor(kGray+1);
  pipe -> SetFillColor(kGray+1);
  pipe -> SetShellDensity(GetDensity("Polyvinyltoluene"));
  pipe -> SetHollowDensity(GetDensity("Water"));
  pipe -> SetAtomicNumber(GetAtomicNumber("Polyvinyltoluene"));
  pipe -> SetAtomicMass(GetAtomicMass("Polyvinyltoluene"));
  pipe -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("Polyvinyltoluene"));
  pipe -> SetSternheimerParameters(GetSternheimerParameters("Polyvinyltoluene"));

  RectangularPyramid* anomaly = new RectangularPyramid(0.0, 0.0, 185.0, 200.0, 20.0, 30.0);
  anomaly -> RotateTheta(TMath::Pi());
  anomaly -> MoveXYZ(100.0, 0.0, 0.0);
  anomaly -> SetLineColor(kAzure-4);
  anomaly -> SetFillColor(kAzure-4);
  anomaly -> SetDensity(GetDensity("Water"));
  anomaly -> SetAtomicNumber(GetAtomicNumber("Water"));
  anomaly -> SetAtomicMass(GetAtomicMass("Water"));
  anomaly -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("Water"));
  anomaly -> SetSternheimerParameters(GetSternheimerParameters("Water"));
//  anomaly -> SetDensity(GetDensity("Air"));

  Pipe* manhole0 = new Pipe(detector_position.X()+500.0,detector_position.Y()-75.0, 260.0, 120.0, 60.0, 65.0);
  manhole0 -> SetLineColor(kGray+1);
  manhole0 -> SetFillColor(kGray);
  manhole0 -> SetShellDensity(GetDensity("ShieldingConcrete"));
  manhole0 -> SetHollowDensity(GetDensity("Air"));
  manhole0 -> SetAtomicNumber(GetAtomicNumber("ShieldingConcrete"));
  manhole0 -> SetAtomicMass(GetAtomicMass("ShieldingConcrete"));
  manhole0 -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("ShieldingConcrete"));
  manhole0 -> SetSternheimerParameters(GetSternheimerParameters("ShieldingConcrete"));

  Cylinder* manhole1 = new Cylinder(detector_position.X()+500.0,detector_position.Y()-75.0, 319, 2.0, 60.0);
  manhole1 -> SetLineColor(kGray+3);
  manhole1 -> SetFillColor(kGray+2);
  manhole1 -> SetDensity(GetDensity("Iron"));
  manhole1 -> SetAtomicNumber(GetAtomicNumber("Iron"));
  manhole1 -> SetAtomicMass(GetAtomicMass("Iron"));
  manhole1 -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("Iron"));
  manhole1 -> SetSternheimerParameters(GetSternheimerParameters("Iron"));

  RectangularPrism* manhole2 = new RectangularPrism(0.0,detector_position.Y()-75.0, 100.0, geometry_rangexy, 200.0, 200.0);
  manhole2 -> SetLineColor(kGray+1);
  manhole2 -> SetFillColor(kGray);
  manhole2 -> SetDensity(GetDensity("Air"));

  SciBar50x50* detector = new SciBar50x50(0.0,0.0,0.0);
  detector -> MoveXYZ(detector_position.X(),detector_position.Y(),detector_position.Z());
//  detector -> RotateTheta(TMath::Pi()/3.0);
//  detector -> RotatePhi(TMath::Pi()/4.0);
  detector -> SetLineColor(3);
  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  ground   -> SetPriority(10);
  road     -> SetPriority(9);
  road2    -> SetPriority(8);
  house    -> SetPriority(7);
  manhole2 -> SetPriority(6);
  manhole0 -> SetPriority(5);
  manhole1 -> SetPriority(4);
  pipe     -> SetPriority(3);
  anomaly  -> SetPriority(2);
  detector -> SetPriority(1);
  std::vector<Shape*> shapes;

  shapes.push_back(ground);
  shapes.push_back(road);
  shapes.push_back(road2);
  shapes.push_back(manhole2);
  shapes.push_back(manhole0);
  shapes.push_back(manhole1);
  shapes.push_back(pipe);
  shapes.push_back(anomaly);
  shapes.push_back(detector);
  shapes.push_back(house);

  // 4. Setup Muon Generator
  ReadExpacs* reader = new ReadExpacs(threshold);

  Muon* muon = new Muon(reader);
  muon->SetReferenceZ(detector_position.Z()-15.0);
  muon->SetEndZ(-500.0);
  muon->SetStartZ(2000.0);
  muon->Set1stBoundaryXY(detector_position.X()-25.0, detector_position.X()+25.0,
                         detector_position.Y()-25.0, detector_position.Y()+25.0);

  // 5. Event Loop
  std::vector<Muon*> muon_tracks;
  int n = std::to_string(n_events).length();
  int n_generation_failed = 0;
  int n_below_threshold = 0;
  int i_draw = 0;
  for (int i = 0; i < n_events; ++i) {
    if (i % (n_events / 10) == 0) Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events);

    while (!muon->Generate()) {
      n_generation_failed++;
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

    // Store hits in detector
    std::vector<int> all_hit_ids = detector->IsHitList(*muon);
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
  view->SetRange(-geometry_rangexy/2, -geometry_rangexy/2, 0.0, geometry_rangexy/2, geometry_rangexy/2, 1000.0);
  view->SetView(800.0, 800.0, 0.0, irep);

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
  return SimulationWaterPipe();
}
#endif

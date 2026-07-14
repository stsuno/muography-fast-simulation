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

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

int SimulationKomen(std::string model="ALL") {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

//  std::string output_dir = "data/output/komen_g20pm00p";
  std::string output_dir = "data/output/komen_g00pm20p";
//  std::string output_dir = "data/output/komen_g00pm500p";
//  std::string output_dir = "data/output/komen_nominal";

  double grav =  0.0;
  double mate =  2.0;

  double geometry_rangez  = 2000.0+mate;
  double random_scale = 1.0;

  TVector3 detA_position(0.0,0.0,0.0);
  TVector3 detB_position(2000.0,0.0,0.0+grav);
  TVector3 detC_position(1000.0,0.0,200.0+0.5*grav);
  TVector3 MWPC_position, det1_position, det2_position;

  std::string datafile1, datafile2;

  if (model=="SingleMWPC") { 
    MWPC_position.SetXYZ(2000.0,0.0,2400.0+mate); 
    det1_position = MWPC_position;
    det2_position = MWPC_position;
    datafile1="komen_singleM.dat";
    datafile2="komen_singleM.dat";
  }
  if (model=="SingleA")    { 
    MWPC_position.SetXYZ(2000.0,0.0,2400.0+mate); 
    det1_position = detA_position;
    det2_position = detA_position;
    datafile1="komen_singleA.dat";
    datafile2="komen_singleA.dat";
  }
  if (model=="SingleB")    { 
    MWPC_position.SetXYZ(2000.0,0.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = detB_position;
    datafile1="komen_singleB.dat";
    datafile2="komen_singleB.dat";
  }
  if (model=="SingleC")    { 
    MWPC_position.SetXYZ(2000.0,0.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = detC_position;
    datafile1="komen_singleC.dat";
    datafile2="komen_singleC.dat";
  }

  if (model=="DoubleM1B")  {
    MWPC_position.SetXYZ(2000.0,  0.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM1B_B.dat";
    datafile2="komen_doubleM1B_M.dat";
  }
  if (model=="DoubleM2B")  {
    MWPC_position.SetXYZ(2000.0, 100.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM2B_B.dat";
    datafile2="komen_doubleM2B_M.dat";
  }
  if (model=="DoubleM3B")  {
    MWPC_position.SetXYZ(2000.0,-100.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM3B_B.dat";
    datafile2="komen_doubleM3B_M.dat";
  }
  if (model=="DoubleM4B")  {
    MWPC_position.SetXYZ(2000.0, 200.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM4B_B.dat";
    datafile2="komen_doubleM4B_M.dat";
  }
  if (model=="DoubleM5B")  {
    MWPC_position.SetXYZ(2000.0,-200.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM5B_B.dat";
    datafile2="komen_doubleM5B_M.dat";
  }

  if (model=="DoubleM1C")  {
    MWPC_position.SetXYZ(2000.0,  0.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM1C_C.dat";
    datafile2="komen_doubleM1C_M.dat";
  }
  if (model=="DoubleM2C")  {
    MWPC_position.SetXYZ(2000.0, 100.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM2C_C.dat";
    datafile2="komen_doubleM2C_M.dat";
  }
  if (model=="DoubleM3C")  {
    MWPC_position.SetXYZ(2000.0,-100.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM3C_C.dat";
    datafile2="komen_doubleM3C_M.dat";
  }
  if (model=="DoubleM4C")  {
    MWPC_position.SetXYZ(2000.0, 200.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM4C_C.dat";
    datafile2="komen_doubleM4C_M.dat";
  }
  if (model=="DoubleM5C")  {
    MWPC_position.SetXYZ(2000.0,-200.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM5C_C.dat";
    datafile2="komen_doubleM5C_M.dat";
  }

  if (model=="DoubleM1A")  {
    MWPC_position.SetXYZ(2000.0,  0.0,2400.0+mate); 
    det1_position = detA_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM1A_A.dat";
    datafile2="komen_doubleM1A_M.dat";
  }
  if (model=="DoubleM2A")  {
    MWPC_position.SetXYZ(2000.0, 100.0,2400.0+mate); 
    det1_position = detA_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM2A_A.dat";
    datafile2="komen_doubleM2A_M.dat";
  }
  if (model=="DoubleM3A")  {
    MWPC_position.SetXYZ(2000.0,-100.0,2400.0+mate); 
    det1_position = detA_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM3A_A.dat";
    datafile2="komen_doubleM3A_M.dat";
  }
  if (model=="DoubleM4A")  {
    MWPC_position.SetXYZ(2000.0, 200.0,2400.0+mate); 
    det1_position = detA_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM4A_A.dat";
    datafile2="komen_doubleM4A_M.dat";
  }
  if (model=="DoubleM5A")  {
    MWPC_position.SetXYZ(2000.0,-200.0,2400.0+mate); 
    det1_position = detA_position;
    det2_position = MWPC_position;
    datafile1="komen_doubleM5A_A.dat";
    datafile2="komen_doubleM5A_M.dat";
  }

  if (model=="DoubleBC")   { 
    MWPC_position.SetXYZ(2000.0,  0.0,2400.0+mate); 
    det1_position = detB_position;
    det2_position = detC_position;
    datafile1="komen_doubleBC_B.dat";
    datafile2="komen_doubleBC_C.dat";
  }
  if (model=="DoubleCA")   { 
    MWPC_position.SetXYZ(2000.0,  0.0,2400.0+mate); 
    det1_position = detC_position;
    det2_position = detA_position;
    datafile1="komen_doubleCA_C.dat";
    datafile2="komen_doubleCA_A.dat";
  }

  std::string pdffile      = "komen.pdf";

  double threshold = 1000.0;
  double flux = 0.008235+0.007109;     // muon threshold=0
  if (threshold==1000.0) { flux = 0.006446+0.005470; }

  double thr_muon = 1000.0;
  if (model=="SingleB"   || model=="SingleC"   
   || model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B" 
   || model=="DoubleM4B" || model=="DoubleM5B"
   || model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C" 
   || model=="DoubleM4C" || model=="DoubleM5C") {
    thr_muon=4000.0; 
    flux = 0.003288+0.002752;
  }
  if (model=="DoubleCA") {
    thr_muon=50000.0; 
    flux = 0.000126+0.000106;
  }

  const double Mmu = 105.658;    // muon mass [MeV]

  int n_events = flux*100.0*100.0*60*60*24;  // 1-day events

  if (model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B"
   || model=="DoubleM4B" || model=="DoubleM5B") {
    n_events = flux*100.0*100.0*60*60*24*228/1511049;  // 787 events per day
  }
  if (model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
   || model=="DoubleM4C" || model=="DoubleM5C") {
    n_events = flux*100.0*100.0*60*60*24*258/1877405;  // 717 events per day
  }
  if (model=="DoubleM1A" || model=="DoubleM2A" || model=="DoubleM3A"
   || model=="DoubleM4A" || model=="DoubleM5A") {
    n_events = flux*100.0*100.0*60*60*24*114/4390472;  // 344 events per day
  }
  if (model=="DoubleBC") {
    n_events = flux*100.0*100.0*60*60*24*869/24689698;  // 466 events per day
  }
  if (model=="DoubleCA") {
    n_events = flux*100.0*100.0*60*60*24*1/1000;       // 200 events per day
  }

  n_events *= 90.0;     // 90 days

  double weight = 1.0;
  if (model=="SingleMWPC") { weight = 1.0/100.0/90.0; }    // number of events per day : 10295423
  if (model=="SingleA")    { weight = 1.0/100.0/90.0; }    //                          :  5218560
  if (model=="SingleB")    { weight = 1.0/50.0/90.0; }     //                          :  5218560
  if (model=="SingleC")    { weight = 1.0/50.0/90.0; }     //                          :  5218560
  n_events *= weight;

  int n_draw = 20;
  TRandom *gaus = new TRandom(231517);    // multiple scattering

  // 2. Setup Output
  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    ::Info("Simulation", "Created output directory: %s", output_dir.c_str());
  }

  // 3. Setup Geometry
  // O:47%, Si:28%, Al:8%, Fe:5%, Ca:4%, Na:3%, K:3%, Mg:2%
  const double atomic_soil  = 0.47* 8+0.28*14+0.08*13+0.05*26+0.04*20+0.03*11+0.03*19+0.02*12;
  const double mass_soil    = 0.47*16+0.28*28+0.08*27+0.05*56+0.04*40+0.03*23+0.03*39+0.02*24;
  const double density_soil = GetDensity("NormalSoil");
  const double SoilX0 = 716.4*mass_soil/(atomic_soil*(atomic_soil+1.0)*TMath::Log(287.0/TMath::Sqrt(atomic_soil)))/density_soil;   // Radiation length

  RectangularPrism* ground = new RectangularPrism(4000.0, 0.0, 0.5*geometry_rangez+300.0, 
                                                  5000.0, 50000.0, geometry_rangez);
  ground -> SetLineColor(791);
  ground -> SetFillColor(791);
  ground -> SetDensity(density_soil);
  ground -> SetAtomicNumber(GetAtomicNumber("StandardRock"));
  ground -> SetAtomicMass(GetAtomicMass("StandardRock"));
  ground -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("StandardRock"));
  ground -> SetSternheimerParameters(GetSternheimerParameters("StandardRock"));

  RectangularPrism* ground2 = new RectangularPrism(1000.0, 0.0, 0.25*geometry_rangez+300.0, 
                                                   1000.0, 50000.0, 0.5*geometry_rangez);
  ground2 -> SetLineColor(791);
  ground2 -> SetFillColor(791);
  ground2 -> SetDensity(density_soil);
  ground2 -> SetAtomicNumber(GetAtomicNumber("StandardRock"));
  ground2 -> SetAtomicMass(GetAtomicMass("StandardRock"));
  ground2 -> SetMeanExcitationEnergy(GetMeanExcitationEnergy("StandardRock"));
  ground2 -> SetSternheimerParameters(GetSternheimerParameters("StandardRock"));

  SciBar100x100* detA = new SciBar100x100(0.0,0.0,0.0);
  detA -> MoveXYZ(detA_position.X(),detA_position.Y(),detA_position.Z());
  detA -> RotateTheta(TMath::Pi()/2.0);
  detA -> RotatePhi(TMath::Pi()/2.0);
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

  SciBar100x100* detC = new SciBar100x100(0.0,0.0,0.0);
  detC -> MoveXYZ(detC_position.X(),detC_position.Y(),detC_position.Z());
  detC -> SetLineColor(3);
  detC -> SetFillColor(3);
  detC -> SetDensity(0.0);

  SciBar100x100* MWPC = new SciBar100x100(0.0,0.0,0.0);
  MWPC -> MoveXYZ(MWPC_position.X(),MWPC_position.Y(),MWPC_position.Z());
  MWPC -> SetLineColor(3);
  MWPC -> SetFillColor(3);
  MWPC -> SetDensity(0.0);

  MWPC    -> SetPriority(1);
  detB    -> SetPriority(2);
  detC    -> SetPriority(3);
  detA    -> SetPriority(4);
  ground  -> SetPriority(5);
  ground2 -> SetPriority(6);
  std::vector<Shape*> shapes;

  shapes.push_back(ground2);
  shapes.push_back(ground);
  shapes.push_back(MWPC);
  shapes.push_back(detB);
  shapes.push_back(detC);
  shapes.push_back(detA);

  // 4. Setup Muon Generator
//  ReadExpacs* reader = new ReadExpacs(threshold);
  ReadExpacs* reader = new ReadExpacs(thr_muon);

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
  std::vector<Muon*> muon_tracks;

  //==============================================================================
  // Simulate detector single A,B,MWPC
  //==============================================================================
  std::string output1_path = output_dir + "/" + datafile1;
  std::ofstream output1_dat(output1_path);
  if (!output1_dat.is_open()) {
    SysError("Simulation::Initialize", "Could not open output file: %s", datafile1.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  std::string output2_path = output_dir + "/" + datafile2;
  std::ofstream output2_dat(output2_path);
  if (!output2_dat.is_open()) {
    SysError("Simulation::Initialize", "Could not open output file: %s", datafile2.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }

  Muon* muon = new Muon(reader, 123456*random_scale);
  muon->SetReferenceZ(det1_position.Z());
  muon->SetEndZ(-500.0);
  muon->SetStartZ(3000.0);
  muon->Set1stBoundaryXY(det1_position.X()-50.0, det1_position.X()+50.0,
                         det1_position.Y()-50.0, det1_position.Y()+50.0);
  if (model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B"
   || model=="DoubleM4B" || model=="DoubleM5B"
   || model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
   || model=="DoubleM4C" || model=="DoubleM5C"
   || model=="DoubleM1A" || model=="DoubleM2A" || model=="DoubleM3A"
   || model=="DoubleM4A" || model=="DoubleM5A"
   || model=="DoubleBC"  || model=="DoubleCA") {
    muon->SetSecondReferenceZ(det2_position.Z());
    muon->Set2ndBoundaryXY(det2_position.X()-50.0, det2_position.X()+50.0,
                           det2_position.Y()-50.0, det2_position.Y()+50.0);
  }

  // 5. Event Loop
  n = std::to_string(n_events).length();
  n_generation_failed = 0;
  n_below_threshold = 0;
  i_draw = 0;
  for (int i = 0; i < n_events; ++i) {
    if (i%(n_events/10)==0) { 
      Info("Simulation::EventLoop", "Progress: %*d / %d processed", n, i, n_events); 
    }

    while (!muon->Generate()) {
      n_generation_failed++;
    }
    double final_energy = 0.0;

    //=====================================================================
    // Simulate multiple scattering
    //=====================================================================
    if (model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B"
     || model=="DoubleM4B" || model=="DoubleM5B"  
     || model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
     || model=="DoubleM4C" || model=="DoubleM5C") {

      std::vector<TVector3> pathHits = ground->HitPoints(*muon);
      if (pathHits.size()!=2) {
        std::cout << "Strange number of hits in the ground " << pathHits.size() << std::endl;
        abort();
      }

      double x0 = pathHits[0].X();
      double y0 = pathHits[0].Y();
      double z0 = pathHits[0].Z();

      double x1 = pathHits[1].X();
      double y1 = pathHits[1].Y();
      double z1 = pathHits[1].Z();

      if (model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
       || model=="DoubleM4C" || model=="DoubleM5C") {
        std::vector<TVector3> pathHits2 = ground2->HitPoints(*muon);
        if (pathHits2.size()!=2) {
          std::cout << "Strange number of hits in the ground2 " << pathHits2.size() << std::endl;
          abort();
        }
        x1 = pathHits2[1].X();
        y1 = pathHits2[1].Y();
        z1 = pathHits2[1].Z();
      }

      double distL = TMath::Sqrt(TMath::Power(x0-x1,2)+TMath::Power(y0-y1,2)+TMath::Power(z0-z1,2));

      double Emu = muon->GetEnergyMuon();    // [MeV]
      double pc   = std::sqrt(Emu*Emu-Mmu*Mmu);
      double beta = pc/Emu;
      double scatter_angle = 13.6/beta/pc*std::sqrt(distL/SoilX0)*(1.0+0.038*std::log(distL/SoilX0));
      double smeared_angle = gaus->Gaus(0.0,scatter_angle);
      double smeared_phi   = 2.0*TMath::Pi()*gaus->Rndm();

      double theta = muon->GetTheta()+smeared_angle;
      double phi   = muon->GetPhi();
      double depth = TMath::Abs(z0-z1);
      double tanTh = TMath::Tan(theta);

      double x1new = x0-depth*tanTh*TMath::Cos(phi);
      double y1new = y0-depth*tanTh*TMath::Sin(phi);
      double z1new = z0-depth;

      double A = (z0-z1new)/(x0-x1new);
      double C = (x0*z1new-x1new*z0)/(x0-x1new);
      double B = (z0-z1new)/(y0-y1new);
      double D = (y0*z1new-y1new*z0)/(y0-y1new);

      double startZ = muon->GetStartZ();
      double endZ   = muon->GetEndZ();

      Muon scattered_muon;
      scattered_muon.SetEnergy(Emu);
      scattered_muon.SetStartPoint(TVector3((startZ-C)/A,(startZ-D)/B,startZ));
      scattered_muon.SetEndPoint(TVector3((endZ-C)/A,(endZ-D)/B,endZ));
      //=====================================================================

      auto [initial_energy, loss_energy] = scattered_muon.CalculateEnergyLoss(shapes);
      final_energy = initial_energy - loss_energy;
    }
    else {
      auto [initial_energy, loss_energy] = muon->CalculateEnergyLoss(shapes);
      final_energy = initial_energy - loss_energy;
    }

    if (final_energy <= threshold) {
      n_below_threshold++;
      continue;
    }

    // Store hits in detector
    std::vector<int> hit1_ids;
    std::vector<int> hit2_ids;

    if (model=="SingleMWPC") { 
      hit1_ids=MWPC->IsHitList(*muon);
      hit2_ids=MWPC->IsHitList(*muon);
    }
    if (model=="SingleA")    { 
      hit1_ids=detA->IsHitList(*muon); 
      hit2_ids=detA->IsHitList(*muon); 
    }
    if (model=="SingleB")    { 
      hit1_ids=detB->IsHitList(*muon); 
      hit2_ids=detB->IsHitList(*muon); 
    }
    if (model=="SingleC")    { 
      hit1_ids=detC->IsHitList(*muon); 
      hit2_ids=detC->IsHitList(*muon); 
    }

    if (model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B"
     || model=="DoubleM4B" || model=="DoubleM5B") {
      hit1_ids=detB->IsHitList(*muon); 
      hit2_ids=MWPC->IsHitList(*muon); 
    }

    if (model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
     || model=="DoubleM4C" || model=="DoubleM5C") {
      hit1_ids=detC->IsHitList(*muon); 
      hit2_ids=MWPC->IsHitList(*muon); 
    }

    if (model=="DoubleM1A" || model=="DoubleM2A" || model=="DoubleM3A"
     || model=="DoubleM4A" || model=="DoubleM5A") {
      hit1_ids=detA->IsHitList(*muon); 
      hit2_ids=MWPC->IsHitList(*muon); 
    }

    if (model=="DoubleBC")   {
      hit1_ids=detB->IsHitList(*muon); 
      hit2_ids=detC->IsHitList(*muon); 
    }
    if (model=="DoubleCA")   {
      hit1_ids=detC->IsHitList(*muon); 
      hit2_ids=detA->IsHitList(*muon); 
    }

    if (hit1_ids.size()<8 || hit2_ids.size()<8) {
      n_nohit_detector++;
      continue;
    }

    for (int id : hit1_ids) {
      output1_dat << std::setw(5) << id;
    }
    output1_dat << "\n";      

    if (model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B"
     || model=="DoubleM4B" || model=="DoubleM5B"
     || model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
     || model=="DoubleM4C" || model=="DoubleM5C"
     || model=="DoubleM1A" || model=="DoubleM2A" || model=="DoubleM3A"
     || model=="DoubleM4A" || model=="DoubleM5A"
     || model=="DoubleBC" || model=="DoubleCA") {
      for (int id : hit2_ids) {
        output2_dat << std::setw(5) << id;
      }
      output2_dat << "\n";      
    }

    if (i_draw < n_draw) {
      Muon* track = new Muon(*muon);
      track -> SetMarkerColor(kPink);
      track -> SetLineColor(kViolet);
      muon_tracks.push_back(track);
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
  output1_dat.close();
  if (model=="DoubleM1B" || model=="DoubleM2B" || model=="DoubleM3B"
   || model=="DoubleM4B" || model=="DoubleM5B"
   || model=="DoubleM1C" || model=="DoubleM2C" || model=="DoubleM3C"
   || model=="DoubleM4C" || model=="DoubleM5C"
   || model=="DoubleM1A" || model=="DoubleM2A" || model=="DoubleM3A"
   || model=="DoubleM4A" || model=="DoubleM5A"
   || model=="DoubleBC" || model=="DoubleCA") {
    output2_dat.close();
  }
  Info("Simulation::Finalize", "Output File: %s", output1_path.c_str());

  // 8. Visualize
  int irep;
  TCanvas* canvas = new TCanvas("canvas", "Simulation", 0, 0, 800, 800);
  TView* view = TView::CreateView();
//  view->SetRange(-2000.0, -2000.0, 0.0, 2000.0, 2000.0, 3000.0);
  view->SetRange(0.0, 0.0, 0.0, 2000.0, 2000.0, 2000.0);
  view->SetView(800.0, 800.0, 0.0, irep);

  for (auto shape : shapes) shape->Draw("f");
  for (auto track : muon_tracks) { track->Draw("m"); }

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

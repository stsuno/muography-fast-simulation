#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "TROOT.h"
#include "TSystem.h"
#include "TAuthenticate.h"
#include "TFile.h"
#include "TTree.h"

#include "KomenAnalysis.h"
#include "SakurajimaAnalysis.h"

int main(int argc, char **argv) {

  int imode = std::stoi(argv[1]);
  std::string FDirectoryName = argv[2];
  std::string inputfilename1 = argv[3];
  std::string inputfilename2;

  //=====================================================
  // imode==1: single detector mode for KomenAnalysis
  //      ==2: double detector mode for KomenAnalysis
  //
  // imode==3: SakurajimaAnalysis
  //=====================================================
  if (imode<1 || imode>3) {
    std::cout << "ERROR::Set mode=1 or 2." << std::endl;
    abort();
  }

  std::cout << "Input filename 1 : " << inputfilename1 << std::endl;
  if (imode==2) {
    inputfilename2=argv[4]; 
    std::cout << "Input filename 2 : " << inputfilename2 << std::endl;
  }

  TChain* tree = new TChain("track_tree");
  TString rin = FDirectoryName + inputfilename1 + "/track_tree";
  tree->AddFile(rin,-1); 

  TString rout = "hist_" + inputfilename1;
  std::cout << "Output File name : " << rout << std::endl;

  if (imode==1) {
    KomenAnalysis myana(tree);
    myana.initialize(rout);
    myana.SetMode(imode);
    myana.execute();
    myana.finalize();
  }
  else if (imode==2) {
    KomenAnalysis myana(tree);
    myana.initialize(rout);
    myana.SetMode(imode);
    myana.SetDirectryName(FDirectoryName);
    myana.SetCheckFile(inputfilename2);

    TVector3 MWPC_position(0.0,0.0,0.0);;
    if (inputfilename1.find("doubleM1B")!=std::string::npos) { MWPC_position.SetXYZ(0.0,   0.0,2400.0); }
    if (inputfilename1.find("doubleM2B")!=std::string::npos) { MWPC_position.SetXYZ(0.0, 100.0,2400.0); }
    if (inputfilename1.find("doubleM3B")!=std::string::npos) { MWPC_position.SetXYZ(0.0,-100.0,2400.0); }
    if (inputfilename1.find("doubleM4B")!=std::string::npos) { MWPC_position.SetXYZ(0.0, 200.0,2400.0); }
    if (inputfilename1.find("doubleM5B")!=std::string::npos) { MWPC_position.SetXYZ(0.0,-200.0,2400.0); }

    if (inputfilename1.find("doubleM1C")!=std::string::npos) { MWPC_position.SetXYZ(1000.0,   0.0,2200.0); }
    if (inputfilename1.find("doubleM2C")!=std::string::npos) { MWPC_position.SetXYZ(1000.0, 100.0,2200.0); }
    if (inputfilename1.find("doubleM3C")!=std::string::npos) { MWPC_position.SetXYZ(1000.0,-100.0,2200.0); }
    if (inputfilename1.find("doubleM4C")!=std::string::npos) { MWPC_position.SetXYZ(1000.0, 200.0,2200.0); }
    if (inputfilename1.find("doubleM5C")!=std::string::npos) { MWPC_position.SetXYZ(1000.0,-200.0,2200.0); }

    if (inputfilename1.find("doubleM1A")!=std::string::npos) { MWPC_position.SetXYZ(2000.0,   0.0,2400.0); }
    if (inputfilename1.find("doubleM2A")!=std::string::npos) { MWPC_position.SetXYZ(2000.0, 100.0,2400.0); }
    if (inputfilename1.find("doubleM3A")!=std::string::npos) { MWPC_position.SetXYZ(2000.0,-100.0,2400.0); }
    if (inputfilename1.find("doubleM4A")!=std::string::npos) { MWPC_position.SetXYZ(2000.0, 200.0,2400.0); }
    if (inputfilename1.find("doubleM5A")!=std::string::npos) { MWPC_position.SetXYZ(2000.0,-200.0,2400.0); }
    myana.SetMWPCPosition(MWPC_position);
    myana.execute();
    myana.finalize();
  }
  else if (imode==3) {
    SakurajimaAnalysis myana(tree);
    myana.initialize(rout);
    myana.execute();
    myana.finalize();
  }
  std::cout << "End of run." << std::endl;

  return 0;
}

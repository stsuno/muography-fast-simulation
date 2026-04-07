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

int main(int argc, char **argv) {
  //================
  // Input & Output
  //================
  TString FDirectoryName;
  char *inputfilename = argv[1];
  std::cout << "Input card name : " << inputfilename << std::endl;

  TString FRootFileName;

  Int_t totFileNumber = -2;
  std::ifstream dataTmp(inputfilename);
  while (!dataTmp.eof()) {
    char stmp[1000];
    dataTmp >> stmp;
    totFileNumber++;
  }

  std::ifstream data(inputfilename);
  data >> FDirectoryName;
  TChain* tree = new TChain("track_tree");

  char s[1000];
  while (!data.eof()) {
    data >> s;
    TString rin = FDirectoryName + s + "/track_tree";
    TString rinA = FDirectoryName + s;
    if (!data.eof()) {
      std::cout << "Input File name  : " << s << std::endl;
      tree->AddFile(rin,-1); 
    }
  }

  FRootFileName = s;
  TString rout = "hist_" + FRootFileName;
  std::cout << "Output File name : " << rout << std::endl;

  KomenAnalysis myana(tree);
  myana.initialize(rout);
  myana.execute();
  myana.finalize();

  std::cout << "End of run." << std::endl;

  return 0;
}

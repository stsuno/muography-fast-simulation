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

#include "ZmumuAnalysis.h"
#include "makeTree.h"
#include "makeEventTree2.h"

using namespace std;

int main(int argc, char **argv) {
  //================
  // Input & Output
  //================
  Int_t iFileNumber = 0;

  TString FDirectoryName;
  char *inputfilename = argv[1];
  int irunMode = (int)atoi(argv[2]);
  int iSelection = -1;
  if (irunMode==1) {
    iSelection = (int)atoi(argv[3]);
  }

  cout << "Run mode : " << irunMode << endl;
  cout << "Input card name : " << inputfilename << endl;
  if (irunMode==1) {
    cout << "Selection mode : " << iSelection << endl;
  }

  TString FRootFileName;

  Int_t totFileNumber = -2;
  ifstream dataTmp(inputfilename);
  while (!dataTmp.eof()) {
    char stmp[1000];
    dataTmp >> stmp;
    totFileNumber++;
  }

  ifstream data(inputfilename);
  data >> FDirectoryName;
  TChain* tree = new TChain("eventtree");

  char s[1000];
  while (!data.eof()) {
    data >> s;

    TString rin = FDirectoryName + s + "/eventtree";

    TString rinA = FDirectoryName + s;

    if (!data.eof()) {
      cout << "Input File name  : " << s << endl;
      tree->AddFile(rin,-1); 
    }
  }

  FRootFileName = s;
  TString rout = "hist_" + FRootFileName;
  cout << "Output File name : " << rout << endl;

  if (irunMode==1) { 
    ZmumuAnalysis myana(tree);
    myana.SetSelectionMode(iSelection);
    myana.initialize(rout);
    myana.execute();
    myana.finalize();
  }
  else if (irunMode==2) { 
    makeTree myana(tree);
    myana.initialize(rout);
    myana.execute();
    myana.finalize();
  }
  else if (irunMode==4) { 
    makeEventTree2 myana(tree);
    myana.initialize(rout);
    myana.execute();
    myana.finalize();
  }

  cout << "End of run." << endl;

  return 0;
}

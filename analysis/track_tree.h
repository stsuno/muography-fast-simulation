#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef track_tree_h
#define track_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class track_tree {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t           event_number;
    Int_t           num_hits;
    Double_t        chisq;
    Int_t           ndf;
    Double_t        reduced_chisq;
    Double_t        p_value;
    Double_t        A;
    Double_t        B;
    Double_t        C;
    Double_t        D;
    Double_t        errorA;
    Double_t        errorB;
    Double_t        errorC;
    Double_t        errorD;
    Double_t        phi;
    Double_t        theta;
    std::vector<double>  *residual_x;
    std::vector<double>  *residual_y;

    // List of branches
    TBranch        *b_event_number;   //!
    TBranch        *b_num_hits;   //!
    TBranch        *b_chisq;   //!
    TBranch        *b_ndf;   //!
    TBranch        *b_reduced_chisq;   //!
    TBranch        *b_p_value;   //!
    TBranch        *b_A;   //!
    TBranch        *b_B;   //!
    TBranch        *b_C;   //!
    TBranch        *b_D;   //!
    TBranch        *b_errorA;   //!
    TBranch        *b_errorB;   //!
    TBranch        *b_errorC;   //!
    TBranch        *b_errorD;   //!
    TBranch        *b_phi;   //!
    TBranch        *b_theta;   //!
    TBranch        *b_residual_x;   //!
    TBranch        *b_residual_y;   //!

    track_tree(TTree *tree=0);
    virtual ~track_tree();
    virtual void     Init(TTree *tree);
    virtual bool     Notify();

protected:
    virtual void initialize() {};
    virtual void execute() {};
    virtual void finalize() {};
};

#endif

#ifdef track_tree_cxx
track_tree::track_tree(TTree *tree) : fChain(0) {
  Init(tree);
}

track_tree::~track_tree() {
}

void track_tree::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  residual_x = 0;
  residual_y = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
  fChain->SetBranchAddress("num_hits", &num_hits, &b_num_hits);
  fChain->SetBranchAddress("chisq", &chisq, &b_chisq);
  fChain->SetBranchAddress("ndf", &ndf, &b_ndf);
  fChain->SetBranchAddress("reduced_chisq", &reduced_chisq, &b_reduced_chisq);
  fChain->SetBranchAddress("p_value", &p_value, &b_p_value);
  fChain->SetBranchAddress("A", &A, &b_A);
  fChain->SetBranchAddress("B", &B, &b_B);
  fChain->SetBranchAddress("C", &C, &b_C);
  fChain->SetBranchAddress("D", &D, &b_D);
  fChain->SetBranchAddress("errorA", &errorA, &b_errorA);
  fChain->SetBranchAddress("errorB", &errorB, &b_errorB);
  fChain->SetBranchAddress("errorC", &errorC, &b_errorC);
  fChain->SetBranchAddress("errorD", &errorD, &b_errorD);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
  fChain->SetBranchAddress("theta", &theta, &b_theta);
  fChain->SetBranchAddress("residual_x", &residual_x, &b_residual_x);
  fChain->SetBranchAddress("residual_y", &residual_y, &b_residual_y);
  Notify();
}

bool track_tree::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return true;
}

#endif // #ifdef track_tree_cxx

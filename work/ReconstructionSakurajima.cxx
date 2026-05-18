#ifdef __CLING__
R__ADD_INCLUDE_PATH(.);
#endif

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <bitset>

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TEnv.h>
#include <TMath.h>
#include <TError.h>
#include <TSystem.h>

#include "Simulation/SakurajimaMWPC.h"

#include "Reconstruction/ClusterMaker.h"
#include "Reconstruction/LineFitter.h"
#include "Reconstruction/TrackFinder.h"
#include "Reconstruction/Track.h"

#ifndef G_DEBUG_LEVEL
#define G_DEBUG_LEVEL
int g_debug_level = 0;
int g_log_level   = kInfo;
#endif

bool inRange(const TString& timestamp, const TString& start, const TString& end) {
    return (timestamp >= start && timestamp <= end);
}

void parseLogLine(const string& line, vector<vector<int>>& xRounds, vector<vector<int>>& yRounds,  TString &timestampOut) {   
  std::istringstream iss(line);
  std::string token;

  std::string index, T, timestamp;
  int someId;

  std::string AdcLabel;
  std::vector<int> adcValues;

  std::string ThpLabel;
  double thp1 = 0, thp2 = 0, thp3 = 0;
  int thp4 = 0;

  std::string CountLabel;
  int count1 = 0, count2 = 0;
  std::string CheckLabel;
  std::vector<int> checkValues;

  std::string PLabel;
  int pValue = 0;
  std::string HvLabel;
  int hv1 = 0, hv2 = 0;
  double hvFloat = 0.0;

  std::string ChLabel;
  std::vector<std::vector<int>> groups;

  // ---- 先頭部分 ----
  iss >> index;
  iss >> T >> timestamp >> someId;

  // ★ timestamp を呼び出し元へ返す
  timestampOut = timestamp;

  // ---- Adc 部分 ----
  iss >> AdcLabel;  
  while (iss >> token) {
    if (token == "Thp") {
      ThpLabel = token;
      break;
    }
    try {
      adcValues.push_back(std::stoi(token));
    } catch (...) {}
  }

  // ---- Thp 部分 ----
  iss >> thp1 >> thp2 >> thp3 >> thp4;

  // ---- Count 部分 ----
  iss >> CountLabel >> count1 >> count2;

  // ---- Check 部分 ----
  iss >> CheckLabel;
  while (iss >> token) {
    if (token == "P") {
      PLabel = token;
      break;
    }
    try {
      checkValues.push_back(std::stoi(token));
    } catch (...) {}
  }

  // ---- P / Hv 部分 ----
  iss >> pValue >> HvLabel >> hv1 >> hv2 >> hvFloat;

  // ---- Ch 部分 ----
  iss >> ChLabel;

  while (iss >> token) {
    try {
      int count = std::stoi(token);
      std::vector<int> values;
      for (int i = 0; i < count; ++i) {
        int val;
        if (iss >> val) values.push_back(val);
      }
      groups.push_back(values);
    } catch (...) {
      break;
    }
  }

  // ---- 結果を格納 ----
  for (int i = 0, round = 0; i + 1 < groups.size() && round < 8; i += 2, ++round) {
    xRounds[round] = groups[i];
    yRounds[round] = groups[i + 1];
  }

}


int ReconstructionSakurajima() {
  // 1. Load Configuration
  // Apply Log Level Setting
  gErrorIgnoreLevel = g_log_level;

  std::string output_dir = "data/output/sakurajima";

  std::string datafile = "sakurajima.dat";
  std::string rootfile = "sakurajima.root";
  std::string pdffile  = "sakurajima.pdf";

  TVector3 detector_position(0.0,0.0,0.0);

  // 2. Setup Output
  std::string output_path = output_dir + "/" + rootfile;
  if (gSystem->AccessPathName(output_dir.c_str())) {
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    Info("Reconstruction", "Created output directory: %s", output_dir.c_str());
  }
  TFile* output_root = new TFile(output_path.c_str(), "RECREATE");

  // 3. Define TTree
  TTree* track_tree = new TTree("track_tree", "Reconstructed Tracks");
  double m_A, m_B, m_C, m_D;
  double m_Aerr, m_Berr, m_Cerr, m_Derr;
  double m_chisq, m_phi, m_theta;
  int m_num_hits;
  std::vector<double>* m_residual_x = new std::vector<double>();
  std::vector<double>* m_residual_y = new std::vector<double>();
  track_tree->Branch("num_hits",      &m_num_hits,      "num_hits/I");
  track_tree->Branch("chisq",         &m_chisq,         "chisq/D");
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

  // 5. Initialize Reconstruction Tools
  SakurajimaMWPC* detector = new SakurajimaMWPC();
  detector -> SetMMOS(8);
  detector -> BuildMWPC(0.0,0.0,0.0);
  detector -> MoveXYZ(detector_position.X(),detector_position.Y(),detector_position.Z());
  detector -> SetLineColor(3);
  detector -> SetFillColor(3);
  detector -> SetDensity(0.0);

  ClusterMaker<SakurajimaMWPC*> cluster_maker(detector);
  cluster_maker.SetClusterModel(ClusterMaker<SakurajimaMWPC*>::CombinedLocal);

  TrackFinder track_finder(detector->GetNumLayers(), 1000000.0);
  track_finder.SetMinimumHits(3);

  std::map<std::string, double> start_values;
  start_values["A"] = 0.0;
  start_values["B"] = 0.0;
  start_values["C"] = detector_position.X();
  start_values["D"] = detector_position.Y();
  track_finder.SetStartValues(start_values);

  TString startTime, endTime;
  cout << "開始日時(YYYYMMDDhhmmss)を入力: "; cin >> startTime;
  cout << "終了日時(YYYYMMDDhhmmss)を入力: "; cin >> endTime;

  // ---- 対応パターン ----
  //vector<TString> patterns = {"Mmos01_DT"};
  vector<TString> patterns = {"Mmos08_DT"};
  //vector<TString> patterns = {"Mmos01_DT", "Mmos02_DT"};
  string pattern_numder = "Mmos08_DT";
  // ---- ディレクトリ内の対象ファイルを取得 ----
  vector<TString> files;
  void* dir = gSystem->OpenDirectory(".");
  const char* entry;
  while ((entry = gSystem->GetDirEntry(dir))) {
    TString fname = entry;
    for (auto& pat : patterns) {
      if (fname.BeginsWith(pat) && fname.EndsWith(".ebe")) {
        // タイムスタンプ抽出
        int dtPos = fname.Index("DT");
        if (dtPos >= 0 && fname.Length() >= dtPos + 16) {
          TString timestamp = fname(dtPos + 2, 14);
          if (inRange(timestamp, startTime, endTime)) {
            files.push_back(fname);
          }
        }
        break;
      }
    }
  }
  gSystem->FreeDirectory(dir);

  // ---- 日付順にソート ----
  sort(files.begin(), files.end(), [](const TString &a, const TString &b){
      TString ta = a(a.Index("DT")+2, 14);
      TString tb = b(b.Index("DT")+2, 14);
      return ta < tb;
      });

  // ---- 出力準備 ----
  TString prevDate = "";
  ofstream outFile;
  ofstream countOut("eventcount.txt", ios::trunc);
  if (!countOut.is_open()) {
    cerr << "Error: eventcount.txt を開けませんでした\n";
    return 0;
  }
  int c0=0, c1=0, c2=0, c3=0, c4=0, c5=0, c6=0, c7=0, c8=0;

  // ---- ファイル順に処理 ----
  for (auto& fname : files) {
    cout << "Processing file: " << fname << endl;

    int dtPos = fname.Index("DT");
    TString timestamp = fname(dtPos + 2, 14);
    TString fileDate = timestamp(0,8); // YYYYMMDD

    if (fileDate != prevDate) {
      if (prevDate.Length() != 0) {
        countOut << string(prevDate.Data()) << " "
                 << c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "
                 << c5<<" "<<c6<<" "<<c7<<" "<<c8 << "\n";
        countOut.flush();
      }
      if (outFile.is_open()) {
        outFile.close();
      }
      string resultFileName = "results_" + pattern_numder + "_" + string(fileDate.Data()) + ".txt";
      outFile.open(resultFileName, ios::trunc);
      if (!outFile.is_open()) {
        cerr << "Error: ファイルを開けませんでした " << resultFileName << endl;
        continue;
      }
      prevDate = fileDate;
      c0=c1=c2=c3=c4=c5=c6=c7=c8=0;
    }

    // ---- ファイル読み込み ----
    ifstream infile(fname);
    if (!infile) {
      cerr << "Error: ファイルを開けませんでした " << fname << endl;
      continue;
    }

    std::string line;
    std::vector<vector<int>> xRounds(8), yRounds(8);
    while (getline(infile, line)) {
      // 初期化
      for (int i = 0; i < 8; i++) { 
        xRounds[i].clear(); yRounds[i].clear(); 
      }

      TString lineTimestamp;
      parseLogLine(line,xRounds,yRounds,lineTimestamp);
      TString lineData = lineTimestamp(0,4) + lineTimestamp(5,2) + lineTimestamp(8,2);

      // ファイル変更時の処理
      if (lineData != prevDate) {
        if (prevDate.Length() != 0) {
          countOut << string(prevDate.Data()) << " "
            << c0<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "
            << c5<<" "<<c6<<" "<<c7<<" "<<c8 << "\n";
          countOut.flush();
        }
        if (outFile.is_open()) {
          outFile.close();
        }
        string resultFileName = "results_" + pattern_numder + "_" + string(lineData.Data()) + ".txt";
        outFile.open(resultFileName, ios::trunc);
        if (!outFile.is_open()) { 
          cerr << "Error: ファイルを開けませんでした" << resultFileName << endl; continue; 
        }
        prevDate = lineData;
        c0=c1=c2=c3=c4=c5=c6=c7=c8=0;
      }

      // 空ベクタチェック
      int emptyCount = 0;
      for (int i = 0; i < 8; i++) {
        if (xRounds[i].empty() || yRounds[i].empty()) emptyCount++;
      }
      switch(emptyCount) {
        case 0: ++c0; break; case 1: ++c1; break; case 2: ++c2; break; 
        case 3: ++c3; break; case 4: ++c4; break; case 5: ++c5; break; 
        case 6: ++c6; break; case 7: ++c7; break; case 8: ++c8; break;
      }
      if (emptyCount > 2) {
        continue;
      }

      int sumhits = 0;
      for (int i=0; i<8; i++) { sumhits+=xRounds[i].size()+yRounds[i].size(); }
      if (sumhits>100) { continue; }



      std::vector<int> hit_ids;
      for (size_t i=0; i<(size_t)xRounds[0].size(); i++) { hit_ids.push_back(xRounds[0][i]+64*0); }
      for (size_t i=0; i<(size_t)yRounds[0].size(); i++) { hit_ids.push_back(yRounds[0][i]+64+1); }
      for (size_t i=0; i<(size_t)xRounds[1].size(); i++) { hit_ids.push_back(xRounds[1][i]+64*2); }
      for (size_t i=0; i<(size_t)yRounds[1].size(); i++) { hit_ids.push_back(yRounds[1][i]+64*3); }
      for (size_t i=0; i<(size_t)xRounds[2].size(); i++) { hit_ids.push_back(xRounds[2][i]+64*4); }
      for (size_t i=0; i<(size_t)yRounds[2].size(); i++) { hit_ids.push_back(yRounds[2][i]+64*5); }
      for (size_t i=0; i<(size_t)xRounds[3].size(); i++) { hit_ids.push_back(xRounds[3][i]+64*6); }
      for (size_t i=0; i<(size_t)yRounds[3].size(); i++) { hit_ids.push_back(yRounds[3][i]+64*7); }
      for (size_t i=0; i<(size_t)xRounds[4].size(); i++) { hit_ids.push_back(xRounds[4][i]+64*8); }
      for (size_t i=0; i<(size_t)yRounds[4].size(); i++) { hit_ids.push_back(yRounds[4][i]+64*9); }
      for (size_t i=0; i<(size_t)xRounds[5].size(); i++) { hit_ids.push_back(xRounds[5][i]+64*10); }
      for (size_t i=0; i<(size_t)yRounds[5].size(); i++) { hit_ids.push_back(yRounds[5][i]+64*11); }
      for (size_t i=0; i<(size_t)xRounds[6].size(); i++) { hit_ids.push_back(xRounds[6][i]+64*12); }
      for (size_t i=0; i<(size_t)yRounds[6].size(); i++) { hit_ids.push_back(yRounds[6][i]+64*13); }
      for (size_t i=0; i<(size_t)xRounds[7].size(); i++) { hit_ids.push_back(xRounds[7][i]+64*14); }
      for (size_t i=0; i<(size_t)yRounds[7].size(); i++) { hit_ids.push_back(yRounds[7][i]+64*15); }

      // Step A: Clustering (Raw hits -> SpacePoints)
      std::vector<SpacePoint> all_space_points = cluster_maker.Execute(hit_ids);

      std::vector<SpacePoint> space_points;
      for (const auto& space_point : all_space_points) {
        SpacePoint sp(space_point.X(), space_point.Y(), space_point.Z());
        sp.SetErrorX(space_point.GetErrorX());
        sp.SetErrorY(space_point.GetErrorY());
        sp.SetErrorZ(space_point.GetErrorZ());
        sp.SetUnitIndex(space_point.GetUnitIndex());
        space_points.push_back(sp);
      }

      // Step B: Track Finding (Combinatorial Search and Fitting)
      track_finder.GetLineFit().SetFittingParameters(&start_values, nullptr, nullptr);
      Track best_track = track_finder.GetBestTrack(space_points);

      // Step C: Recording (Fill TTree if a valid track was found)
      if (best_track.IsValid()) {
        m_num_hits      = best_track.GetNumHits();
        m_chisq         = best_track.GetChiSquare();
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
  }
  output_root->Write();
  output_root->Close();
  Info("Reconstruction", "Saved histograms to %s", output_path.c_str());

  return 0;
}


#ifndef __CLING__
int main(int argc, char** argv) {
  return ReconstructionSakurajima();
}
#endif

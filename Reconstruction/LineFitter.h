#ifndef LINE_FITTER_H
#define LINE_FITTER_H

#include <map>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include <TMinuit.h>
#include <TMath.h>
#include <TError.h>

#include "Reconstruction/SpacePoint.h"

// Structure to hold fit status and results
struct FitResult {
  bool flag_migrad = false;
  bool flag_hesse = false;
  double chisq = -1.0;
  double reduced_chisq = -1.0;
  double p_value = 0.0;
  int ndf = 0;
  std::map<std::string, double> values;
  std::map<std::string, double> errors;
};

class LineFitter : public TObject {
  public:
    LineFitter(const std::string& mode = "vertical");
    virtual ~LineFitter() {}

    // Minimize chisquare using MIGRAD and HESSE.
    FitResult ExecuteFit();

    // Add space point to the fitter
    void AddSpacePoint(const SpacePoint& space_point) { m_space_points.push_back(space_point); }
    void ClearSpacePoint() { m_space_points.clear(); }

    // Accessors
    void SetFittingParameters(const std::map<std::string, double>* start_values = nullptr,
        const std::map<std::string, double>* step_sizes = nullptr,
        const std::map<std::string, std::pair<double, double>>* value_limits = nullptr);
    std::string GetFittingMode() const { return m_mode; }

    const std::vector<SpacePoint>& GetSpacePoints() const { return m_space_points; }

    // Static Function for TMinuit FCN (Function Chi-Square / Function to minimize)
    static void FCN(int& npar, double* gin, double& f, double* par, int iflag);

  private:
    std::string m_mode = "vertical";
    std::vector<SpacePoint> m_space_points;
    // Fitting configuration
    std::map<std::string, double> m_start_values;
    std::map<std::string, double> m_step_sizes;
    std::map<std::string, std::pair<double, double>> m_value_limits;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline LineFitter::LineFitter(const std::string& mode)
  : m_mode(mode) {
  std::vector<std::string> names = {"A", "B", "C", "D"};
  for (const auto& name : names) {
    bool is_slope = (name == "A" || name == "B");
    m_start_values[name] = 0.0;
    m_step_sizes[name]   = is_slope ? 0.1 : 0.5;
    m_value_limits[name] = std::make_pair(-2000.0, 2000.0);
  }
}

inline FitResult LineFitter::ExecuteFit() {
  FitResult result;
  // for (size_t i = 0; i < m_space_points.size(); ++i) {
  //   Info("LineFitter::ExecuteFit", "Point %lu: (x, y, z) = (%.2f, %.2f, %.2f)", i, m_space_points[i].X(), m_space_points[i].Y(), m_space_points[i].Z());
  //   Info("LineFitter::ExecuteFit", "Error %lu: (x, y, z) = (%.2f, %.2f, %.2f)", i, m_space_points[i].GetErrorX(), m_space_points[i].GetErrorY(), m_space_points[i].GetErrorZ());
  // }
  // Requirement: At least 3 points to have NDF > 0 for 4 parameters

//  if (m_space_points.size() < 3) {
//    ::Error("LineFitter::ExecuteFit", "Not enough space points to fit: %d", static_cast<int>(m_space_points.size()));
//    return {false, false, 0.0};
//  }


  // TMinuit instance (4 parameters: A, B, C, D)
  std::unique_ptr<TMinuit> minuit = std::make_unique<TMinuit>(4);
  minuit->SetFCN(LineFitter::FCN);
  minuit->SetObjectFit(this);

  double arglist[2];
  int ierflg = 0;

  // Set print level (-1: quiet, 0: normal, 1: verbose)
  arglist[0] = -1;
  minuit->mnexcm("SET PRI", arglist, 1, ierflg);

  // Set error definition for Least Squares
//STSTST  arglist[0] = 1.0;
//STSTST  minuit->mnexcm("SET ERR", arglist, 1, ierflg);

  // Define parameters
  const std::vector<std::string> names = {"A", "B", "C", "D"};
  for (int i = 0; i < 4; ++i) {
    const std::string& n = names[i];
    minuit->mnparm(i, n.c_str(), m_start_values[n], m_step_sizes[n], m_value_limits[n].first, m_value_limits[n].second, ierflg);
    // Info("LineFitter::ExecuteFit", "start, step, min, max: %.2f, %.2f, %.2f, %.2f", m_start_values[n], m_step_sizes[n], m_value_limits[n].first, m_value_limits[n].second);
  }

  // // Execute SIMPLEX first to get closer to the minimum
  // arglist[0] = 1000;  // Max calls
  // arglist[1] = 0.1;   // Toletance
  // minuit->mnexcm("SIMPLEX", arglist, 2, ierflg);

  // Execute MIGRAD (minimize fcn) perform minimization
  arglist[0] = 5000;  // Max calls
  arglist[1] = 0.1;   // Toletance
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  bool flag_migrad = (ierflg == 0);
  // Info("LineFitter::ExecuteFit", "ierflg migrad: %d", ierflg);
  // Info("LineFitter::ExecuteFit", "flag migrad: %d", flag_migrad);

  // Execute HESSE (symmetrical error) calculate covariance matrix
  arglist[0] = 5000;  // Max calls
  minuit->mnexcm("HESSE", arglist, 1, ierflg);
  bool flag_hesse = (ierflg == 0);
  // Info("LineFitter::ExecuteFit", "flag hesse: %d", flag_hesse);

  // Retrieve fit statistics
  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  // Store results
  result.flag_migrad = flag_migrad;
  result.flag_hesse = flag_hesse;
  result.chisq = amin;
  // Calculate Reduced Chi-Square and P-Value (NDF = n_points * 2 - n_parameters)
  result.ndf = static_cast<int>(m_space_points.size()) * 2 - 4;
  result.reduced_chisq = (result.ndf > 0) ? amin / result.ndf : -1.0;
  result.p_value = (result.ndf > 0) ? TMath::Prob(amin, result.ndf) : 0.0;

  for (int i = 0; i < 4; ++i) {
    double value, error;
    minuit->GetParameter(i, value, error);
    result.values[names[i]] = value;
    result.errors[names[i]] = error;
  }

  return result;
}

inline void LineFitter::FCN(int& npar, double* gin, double& f, double* par, int iflag) {
  LineFitter* object = static_cast<LineFitter*>(gMinuit->GetObjectFit());
  // TObject* fit_object = gMinuit->GetObjectFit();
  // if (!fit_object) {
  //   ::Error("LineFitter::FCN", "Object fit is null.");
  //   f = 1e10;
  //   return;
  // }
  // LineFitter* object = static_cast<LineFitter*>(fit_object);
  double total_chisq = 0.0;

  bool is_vertical = (object->GetFittingMode() == "vertical");
  bool is_horizontal = (object->GetFittingMode() == "horizontal");
  double A = par[0];
  double B = par[1];
  double C = par[2];
  double D = par[3];

  for (const auto& space_point : object->m_space_points) {
    double obs_x = space_point.X();
    double obs_y = space_point.Y();
    double obs_z = space_point.Z();
    double err_x = space_point.GetErrorX();
    double err_y = space_point.GetErrorY();
    double err_z = space_point.GetErrorZ();

    // Info("LineFitter::FCN", "obs (x, y, z): (%.2f, %.2f, %.2f)", obs_x, obs_y, obs_z);
    // Info("LineFitter::FCN", "err (x, y, z): (%.2f, %.2f, %.2f)", err_x, err_y, err_z);
    if (is_vertical) {
      // Info("LineFitter::FCN", "mode: %s", object->GetFittingMode().c_str());
      double exp_x_z = A * obs_z + C;
      double exp_y_z = B * obs_z + D;
      double delta_x = (obs_x - exp_x_z) / err_x;
      double delta_y = (obs_y - exp_y_z) / err_y;
      total_chisq += delta_x * delta_x + delta_y * delta_y;
//      if (err_x<40.0) {
//        total_chisq += delta_x*delta_x;
//      }
//      if (err_y<40.0) {
//        total_chisq += delta_y*delta_y;
//      }

//      std::cout << exp_x_z << " " << exp_y_z << " " << obs_x << " " << obs_y << " " << obs_z << " " << err_x << " " << err_y << " " << err_z << " " << total_chisq << std::endl;

    }
    else if (is_horizontal) {
      // Info("LineFitter::FCN", "mode: %s", object->GetFittingMode().c_str());
      double exp_z_x = A * obs_x + C;
      double exp_z_y = B * obs_y + D;
      double delta_x = (obs_z - exp_z_x) / err_z;
      double delta_y = (obs_z - exp_z_y) / err_z;
      total_chisq += delta_x * delta_x + delta_y * delta_y;
    }
  }
  f = total_chisq; 
  // Info("LineFitter::FCN", "total chisq: %.2f", total_chisq);
}

inline void LineFitter::SetFittingParameters(const std::map<std::string, double>* start_values,
    const std::map<std::string, double>* step_sizes,
    const std::map<std::string, std::pair<double, double>>* value_limits) {

  if (start_values) {
    for (const auto& it : *start_values) {
      if (m_start_values.count(it.first)) { m_start_values[it.first]=it.second; }
    }
  }
  if (step_sizes) {
    for (const auto& it : *step_sizes) {
      if (m_step_sizes.count(it.first)) { m_step_sizes[it.first]=it.second; }
    }
  }
  if (value_limits) {
    for (const auto& it : *value_limits) {
      if (m_value_limits.count(it.first)) { m_value_limits[it.first]=it.second; }
    }
  }
}

#endif

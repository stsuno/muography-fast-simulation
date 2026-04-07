#ifndef TRACK_FINDER_H
#define TRACK_FINDER_H

#include <map>
#include <string>
#include <vector>

#include <TMath.h>

#include "Reconstruction/SpacePoint.h"
#include "Reconstruction/LineFitter.h"
#include "Reconstruction/Track.h"

class TrackFinder {
  public:
    TrackFinder(int nlayers, double reduced_chisq_threshold);
    virtual ~TrackFinder() {}

    void SetMinimumHits(int minhit) { m_minhit=minhit; }
    void SetStartValues(std::map<std::string, double> start_values) {
      m_line_fitter.SetFittingParameters(&start_values, nullptr, nullptr);
    }

    LineFitter GetLineFit() { return m_line_fitter; }

    // Executes combinatorial search and returns the best track found
    Track GetBestTrack(const std::vector<SpacePoint>& all_space_points);

  private:
    // Internal recursive method with backtracking
    void ExecuteSearch(const std::vector<std::vector<SpacePoint>>& units_points, size_t unit_index, std::vector<SpacePoint>& current_combination);
    int m_num_units, m_minhit;
    LineFitter m_line_fitter;
    double m_reduced_chisq_threshold;

    // Storage for the best candidate found during the search of a single event
    double m_min_reduced_chisq = TMath::Infinity();
    FitResult m_best_result;
    std::vector<SpacePoint> m_best_combination;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline TrackFinder::TrackFinder(int nlayers, double reduced_chisq_threshold)
  : m_num_units(nlayers),
  m_reduced_chisq_threshold(reduced_chisq_threshold) { }

inline Track TrackFinder::GetBestTrack(const std::vector<SpacePoint>& all_space_points) {
  // Reset event-specific variables
  m_min_reduced_chisq = TMath::Infinity();
  m_best_combination.clear();
  m_best_result = FitResult();

  // 1. Group SpacePoints by their unit index
  std::vector<std::vector<SpacePoint>> units_points(m_num_units);
  for (const auto& space_point : all_space_points) {
    int unit_index = space_point.GetUnitIndex();
    if (unit_index >= 0 && unit_index < m_num_units) {
      units_points[unit_index].push_back(space_point);
    }
  }

  // 2. Start recursive search
  std::vector<SpacePoint> current_combination;
  ExecuteSearch(units_points, 0, current_combination);

  // 3. Return the best Track found
  return Track(m_best_result, m_best_combination);
}

inline void TrackFinder::ExecuteSearch(const std::vector<std::vector<SpacePoint>>& units_points, size_t unit_index, std::vector<SpacePoint>& current_combination) {
  // Base case: Reached the last unit
  if (unit_index == units_points.size()) {
    // Requirements: At least 3 points are needed for a meaningful line fit
    if ((int)current_combination.size() < m_minhit) return;

    // Perform fitting for this specific combination
    m_line_fitter.ClearSpacePoint();
    for (const auto& space_point : current_combination) {
      m_line_fitter.AddSpacePoint(space_point);
    }
    FitResult result = m_line_fitter.ExecuteFit();

    if (g_debug_level>0) {
      Info("TrackFitter::ExecuteSearch", "flag (MIGRAD, HESSE): %d, %d", result.flag_migrad, result.flag_hesse);
    }
    // Check if the fit converged and satisfies the chi-square quality cut
    if (result.flag_migrad) {
      if (g_debug_level>0) { Info("TrackFitter::ExecuteSearch", "MIGRAD: ok"); }
      if (result.flag_hesse) {
        if (g_debug_level>0) { Info("TrackFitter::ExecuteSearch", "HESSE: ok"); }
        if (result.reduced_chisq < m_reduced_chisq_threshold) {
          if (g_debug_level>0) { Info("TrackFitter::ExecuteSearch", "Track quality: ok"); }
          if (result.reduced_chisq < m_min_reduced_chisq) {
            m_min_reduced_chisq = result.reduced_chisq;
            m_best_result = result;
            m_best_combination = current_combination;
          }
        } else {
          if (g_debug_level>0) { Info("TrackFitter::ExecuteSearch", "Bad track quality: chisq/ndf = %f", result.reduced_chisq); }
        }
      } else {
        if (g_debug_level>0) { Info("TrackFitter::ExecuteSearch", "HESSE failed to calculate accurate covariance matrix.");	}
      }
    } else {
      if (g_debug_level>0) { Info("TrackFitter::ExecuteSearch", "MIGRAD failed to find minimum."); }
    } 
    return;
  }

  // Branch 1: Try including each SpacePoint from the current unit
  for (const auto& space_point : units_points[unit_index]) {
    current_combination.push_back(space_point);
    ExecuteSearch(units_points, unit_index+1, current_combination);
    current_combination.pop_back();  // Backtrack: remove the point for the next iteration
  }

  // Branch 2: Try skipping the current unit entirely (in case of inefficiency or noise)
  ExecuteSearch(units_points, unit_index+1, current_combination);

}

#endif

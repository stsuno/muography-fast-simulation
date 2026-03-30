#ifndef TRACK_H
#define TRACK_H

#include <limits>
#include <string>
#include <vector>

#include <TMath.h>

#include "SpacePoint.h"
#include "LineFitter.h"

class Track {
  public:
    Track() : m_is_valid(false) {}
    Track(const FitResult& result, const std::vector<SpacePoint>& space_points, const std::string& mode = "vertical");
    virtual ~Track() {}

    // Getter for Hits
    bool IsValid() const { return m_is_valid; }
    const std::vector<SpacePoint>& GetSpacePoints() const { return m_space_points; }
    int GetNumHits() const { return static_cast<int>(m_space_points.size()); }

    // Getter for Fit Statistics
    double GetA() const { return m_a; }
    double GetB() const { return m_b; }
    double GetC() const { return m_c; }
    double GetD() const { return m_d; }
    double GetErrorA() const { return m_aerr; }
    double GetErrorB() const { return m_berr; }
    double GetErrorC() const { return m_cerr; }
    double GetErrorD() const { return m_derr; }
    double GetChiSquare() const { return m_chisq; }
    double GetReducedChiSquare() const { return m_reduced_chisq; }
    double GetPValue() const { return m_p_value; }
    int GetNDF() const { return m_ndf; }

    // Physics Logic
    double GetX(double z) const;
    double GetY(double z) const;
    double GetPhi() const;
    double GetTheta() const;
    std::vector<double> GetResidualX() const;
    std::vector<double> GetResidualY() const;

  private:
    bool m_is_valid = false;
    std::vector<SpacePoint> m_space_points;
    std::string m_mode = "vertical"; // vertical: beam along Z, horizontal: beam along X

    double m_a = TMath::QuietNaN();
    double m_b = TMath::QuietNaN();
    double m_c = TMath::QuietNaN();
    double m_d = TMath::QuietNaN();
    double m_aerr = TMath::QuietNaN();
    double m_berr = TMath::QuietNaN();
    double m_cerr = TMath::QuietNaN();
    double m_derr = TMath::QuietNaN();
    double m_chisq = -1.0;
    double m_reduced_chisq = -1.0;
    double m_p_value = 0.0;
    int m_ndf = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline Track::Track(const FitResult& result, const std::vector<SpacePoint>& space_points, const std::string& mode)
  : m_is_valid(result.flag_migrad), m_space_points(space_points), m_mode(mode) {

  if (m_is_valid) {
    m_a = result.values.at("A");
    m_b = result.values.at("B");
    m_c = result.values.at("C");
    m_d = result.values.at("D");
    m_aerr = result.errors.at("A");
    m_berr = result.errors.at("B");
    m_cerr = result.errors.at("C");
    m_derr = result.errors.at("D");
    m_chisq = result.chisq;
    m_reduced_chisq = result.reduced_chisq;
    m_p_value = result.p_value;
    m_ndf = result.ndf;
  }
}

inline double Track::GetX(double z) const {
  // virtical mode: Linear equation x = A*z + C
  if (m_mode == "vertical") return m_a * z + m_c;
  // horizontal mode: Solve z = A*x + C => x = (z - C) / A
  if (m_mode == "horizontal") return (z - m_c) / m_a;
  return TMath::QuietNaN();
}

inline double Track::GetY(double z) const {
  // virtical mode: Linear equation y = B*z + D
  if (m_mode == "vertical") return m_b * z + m_d;
  // horizontal mode: Solve z = B*y + D => y = (z - D) / B
  if (m_mode == "horizontal") return (z - m_d) / m_b;
  return TMath::QuietNaN();
}

inline double Track::GetPhi() const {
  return TMath::ATan2(m_b, m_a);
}

inline double Track::GetTheta() const {
  double dr_dbeam = TMath::Sqrt(m_a * m_a + m_b * m_b);
  return TMath::ATan(dr_dbeam);
}

inline std::vector<double> Track::GetResidualX() const {
  std::vector<double> residuals;
  for (const auto& space_point : m_space_points) {
    residuals.push_back(space_point.X() - GetX(space_point.Z()));
  }
  return residuals;
}

inline std::vector<double> Track::GetResidualY() const {
  std::vector<double> residuals;
  for (const auto& space_point : m_space_points) {
    residuals.push_back(space_point.Y() - GetY(space_point.Z()));
  }
  return residuals;
}

#endif

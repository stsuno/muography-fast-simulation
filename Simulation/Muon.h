#ifndef MUON_H
#define MUON_H

#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TError.h>

#include "Simulation/ReadExpacs.h"
#include "Simulation/Shape.h"

extern int g_debug_level;

class Muon : public TPolyLine3D {
  public:
    Muon() { };
    Muon(ReadExpacs* reader, UInt_t random_seed = 1);
    virtual ~Muon();

    bool Generate();

    void Set1stBoundaryXY(double xmin, double xmax, double ymin, double ymax);
    void Set2ndBoundaryXY(double xmin, double xmax, double ymin, double ymax);

    std::pair<double, double> CalculateEnergyLoss(const std::vector<Shape*>& shapes);

    // Accessors (Getter and Setter)
    void SetStartZ(double z) { m_start_z = z; }
    double GetStartZ() const { return m_start_z; }
    void SetEndZ(double z) { m_end_z = z; }
    double GetEndZ() const { return m_end_z; }
    void SetReferenceZ(double z) { m_reference_z = z; }
    double GetReferenceZ() const { return m_reference_z; }
    void SetSecondReferenceZ(double z) { m_second_reference_z = z; }
    double GetSecondReferenceZ() const { return m_second_reference_z; }
    void SetCylinderBoundary(double radius, double height, int sign);

    double GetEnergyMuon() const { return m_energy_muon; }
    int    GetCharge() const { return m_charge; }
    double GetTheta()  const { return m_theta; }
    double GetPhi()    const { return m_phi; }

    TVector3 GetStartPoint() const { return TVector3(m_start_x, m_start_y, m_start_z); }
    TVector3 GetEndPoint() const { return TVector3(m_end_x, m_end_y, m_end_z); }
    TVector3 GetReferencePoint() const { return TVector3(m_reference_x, m_reference_y, m_reference_z); }

    void SetEnergyMuon(double energy) { m_energy_muon=energy; }
    void SetStartPoint(TVector3 p) { m_start_x=p.X(); m_start_y=p.Y(); m_start_z=p.Z(); }
    void SetEndPoint(TVector3 p)   { m_end_x  =p.X(); m_end_y  =p.Y(); m_end_z  =p.Z(); }

    // Visualization
    void AddHitMarker(const TVector3& point);
    void Draw(Option_t* option = "") override;
    void SetMarkerColor(Color_t color);

  private:
    // Constants for geometric calculations
    static constexpr double s_epsilon = 1e-9;  // Threshold to determine if a value is effectively zero

    double CalculateSegmentEnergyLoss(double density, double muon_energy, double path_length);

    ReadExpacs* m_reader;
    TRandom3 m_rnd;

    double m_start_x = 0.0, m_start_y = 0.0, m_start_z = 0.0;
    double m_end_x = 0.0, m_end_y = 0.0, m_end_z = 0.0;
    double m_reference_x = 0.0, m_reference_y = 0.0, m_reference_z = 0.0;
    double m_second_reference_z = 0.0;
    double m_theta = 0.0, m_phi = 0.0;
    double m_boundary1_xmin = 0.0, m_boundary1_xmax = 0.0, m_boundary1_ymin = 0.0, m_boundary1_ymax = 0.0;
    double m_boundary2_xmin = 0.0, m_boundary2_xmax = 0.0, m_boundary2_ymin = 0.0, m_boundary2_ymax = 0.0;

    int m_boundary_radsign = 0;
    double m_boundary_radius=0.0, m_boundary_height=0.0;

    double m_energy_muon = 0.0;
    int m_charge = 0;

    TPolyMarker3D* m_hit_markers = nullptr;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline Muon::Muon(ReadExpacs* reader, UInt_t random_seed)
  : TPolyLine3D(2), m_reader(reader), m_rnd(random_seed) {

  m_hit_markers = new TPolyMarker3D();
  m_hit_markers->SetMarkerStyle(8);
  m_hit_markers->SetMarkerSize(0.3);

  this->SetLineWidth(1);
}

inline Muon::~Muon() {
  if (m_hit_markers) delete m_hit_markers;
}

inline void Muon::Set1stBoundaryXY(double xmin, double xmax, double ymin, double ymax) {
  m_boundary1_xmin = xmin, m_boundary1_xmax = xmax, m_boundary1_ymin = ymin, m_boundary1_ymax = ymax;
}

inline void Muon::Set2ndBoundaryXY(double xmin, double xmax, double ymin, double ymax) {
  m_boundary2_xmin = xmin, m_boundary2_xmax = xmax, m_boundary2_ymin = ymin, m_boundary2_ymax = ymax;
}

inline void Muon::SetCylinderBoundary(double radius, double height, int sign) {
  m_boundary_radius  = radius;
  m_boundary_height  = height;
  m_boundary_radsign = sign;
}

inline bool Muon::Generate() {
  // 1. Randomly determine reference point (x, y) on generation plane (1st. boundary)
  if (m_boundary_radius>0.0) {
    double rnd_r   = m_rnd.Rndm();
    double rnd_phi = m_rnd.Rndm();
    double rnd_z   = m_rnd.Rndm();
    m_reference_x = m_boundary_radius*rnd_r*TMath::Cos(2.0*TMath::Pi()*rnd_phi);
    m_reference_y = m_boundary_radius*rnd_r*TMath::Sin(2.0*TMath::Pi()*rnd_phi);
//    m_reference_z = m_boundary_height*(rnd_z-0.5);
  }
  else {
    m_reference_x = m_boundary1_xmin + (m_boundary1_xmax - m_boundary1_xmin) * m_rnd.Rndm();
    m_reference_y = m_boundary1_ymin + (m_boundary1_ymax - m_boundary1_ymin) * m_rnd.Rndm();
  }

  // 2. Determine azimuthal angle (phi = 0 ~ 360 deg.)
  double phi_deg = 360.0 * m_rnd.Rndm();
  double phi = TMath::DegToRad() * phi_deg;
  double cos_phi = TMath::Cos(phi);
  double sin_phi = TMath::Sin(phi);

  // 3. Calculation theta_min/_max for specific phi (Geometric Ray Casting)
  // Distance between 1st. boundary and 2nd. boundary
  double delta_z = TMath::Abs(m_reference_z - m_second_reference_z);
  double theta_min = 0.0;
  double theta_max = 90.0;
  if (m_boundary_radius>0.0) {
    double r_max = m_boundary_radius-TMath::Sqrt(m_reference_x*m_reference_x+m_reference_y*m_reference_y);
    if (m_boundary_radsign>0) {
      theta_min = TMath::Min(TMath::ATan2(r_max,0.5*m_boundary_height-m_reference_z),
                             TMath::ATan2(r_max,0.5*m_boundary_height+m_reference_z))*TMath::RadToDeg();
    }
    else {
      theta_max = TMath::Max(TMath::ATan2(r_max,0.5*m_boundary_height-m_reference_z),
                             TMath::ATan2(r_max,0.5*m_boundary_height+m_reference_z))*TMath::RadToDeg();
    }
  }
  else {
//    if (delta_z > s_epsilon) {
    if (TMath::Abs(m_second_reference_z)>s_epsilon) {
      // Check if reference point is inside or outside 2nd. boundary
      bool is_inside = (m_reference_x >= m_boundary2_xmin && m_reference_x <= m_boundary2_xmax &&
          m_reference_y >= m_boundary2_ymin && m_reference_y <= m_boundary2_ymax);

      // Case1: Reference point is INSIDE
      if (is_inside) {
        // Find nearest boundary in current phi direction
        double x_wall = (cos_phi > s_epsilon) ? m_boundary2_xmax : m_boundary2_xmin;
        double y_wall = (sin_phi > s_epsilon) ? m_boundary2_ymax : m_boundary2_ymin;

        double r_x = (TMath::Abs(cos_phi) > s_epsilon) ? (x_wall - m_reference_x) / cos_phi : 1e+18;
        double r_y = (TMath::Abs(sin_phi) > s_epsilon) ? (y_wall - m_reference_y) / sin_phi : 1e+18;

        // Distance to closest wall
        double r_max = TMath::Min(r_x, r_y);

        // Final calculation of theta_max
        theta_min = 0.0;
        theta_max = TMath::ATan2(r_max, delta_z) * TMath::RadToDeg();
      }
      else {
        // Case2: Reference point is OUTSIDE
        // Perform intersection checks with all 4 edges and find maximum valid distance
        std::vector<double> valid_radius;
        double x_wall[2] = {m_boundary2_xmin, m_boundary2_xmax};
        double y_wall[2] = {m_boundary2_ymin, m_boundary2_ymax};

        // Check intersection with X-boundaries (vertical edges)
        for (double x : x_wall) {
          if (TMath::Abs(cos_phi) > s_epsilon) {
            double radius = (x - m_reference_x) / cos_phi;
            if (radius > s_epsilon) {
              // Check if intersection point is with in Y-range of detector
              double intersection_y = m_reference_y + radius * sin_phi;
              if (intersection_y>=m_boundary2_ymin && intersection_y<=m_boundary2_ymax) {
                valid_radius.push_back(radius);
              }
            }
          }
        }
        // Check intersection with Y-boundaries (horizontal edges)
        for (double y : y_wall) {
          if (TMath::Abs(sin_phi) > s_epsilon) {
            double radius = (y - m_reference_y) / sin_phi;
            if (radius > s_epsilon) {
              // Check if intersection point is with in X-range of detector
              double intersection_x = m_reference_x + radius * cos_phi;
              if (intersection_x>=m_boundary2_xmin && intersection_x<=m_boundary2_xmax) {
                valid_radius.push_back(radius);	    
              }
            }
          }
        }
        // If no walls are intersected, direction is outside field of view
        if (valid_radius.empty()) {
          if (g_debug_level>0) { Info("Muon::Generate", "Valid radius is empty"); }
          return false;
        }

        // Use furthest and nearest intersection points to define the angle range
        auto min_max = std::minmax_element(valid_radius.begin(), valid_radius.end());
        theta_min = TMath::ATan2(*min_max.first, delta_z) * TMath::RadToDeg();
        theta_max = TMath::ATan2(*min_max.second, delta_z) * TMath::RadToDeg();
      }
    }
  }

  // 4. Acceptance-Rejection sampling to eliminate geometric bias
  const std::vector<double>& angles = m_reader->GetSamplingAngles();
  const std::vector<double>& flux_p = m_reader->GetAngularFluxP();
  const std::vector<double>& flux_m = m_reader->GetAngularFluxM();

  std::vector<double> bounded_flux_pm(angles.size(), 0.0);
  double total_visible_flux = 0.0;

  // Summing plus and minus flux first to determine theta before charge
  for (size_t i = 0; i < angles.size(); ++i) {
    if (angles[i] >= theta_min && angles[i] <= theta_max) {
      bounded_flux_pm[i] = flux_p[i] + flux_m[i];
      total_visible_flux += bounded_flux_pm[i];
    }
  }

  // Reference flux: total sky sum (plus + minus) to normalize acceptance probability
  double max_sky_flux = m_reader->GetTotalFluxP() + m_reader->GetTotalFluxM();

  // Rejection: Acceptance probability based on visible 'sky flux' relative to total
  if (total_visible_flux <= s_epsilon || m_rnd.Rndm() > (total_visible_flux /max_sky_flux)) {
    if (g_debug_level>0) { Info("Muon::Generate", "Acceptance rejection"); }
    return false;  // Event rejected to maintain uniform background
  }

  // 5. Theta (zenith angle) sampling via inverse transform method
  int theta_index = -1;
  double angle_probability = 0.0;
  double m_rndangle_probability = m_rnd.Rndm();
  for (size_t i = 0; i < bounded_flux_pm.size(); ++i) {
    // calculate PDF and accumulate to build CDF
    angle_probability += bounded_flux_pm[i] / total_visible_flux;
    if (angle_probability >= m_rndangle_probability) {
      theta_index = static_cast<int>(i);
      break;
    }
  }
  if (theta_index == -1) theta_index = static_cast<int>(bounded_flux_pm.size()) - 1;
  double theta_deg = angles[theta_index]+m_rnd.Rndm();

  m_phi   = TMath::DegToRad()*phi_deg;
  m_theta = TMath::DegToRad()*theta_deg;

  // 6. Charge and energy sampling via inverse transform method with linear interpolation
  // Determine charge based on the plus/minus ratio at the sampled angle
  double prob_plus = flux_p[theta_index] / (flux_p[theta_index] + flux_m[theta_index]);
  m_charge = (m_rnd.Rndm() < prob_plus) ? 1 : -1;

  // Choose appropriate CDF and energy list based on the determined charge
  const auto& target_energy_cdf = (m_charge > 0) ? m_reader->GetEnergyCDFP() : m_reader->GetEnergyCDFM();
  const std::vector<double>& energy_probability = target_energy_cdf[theta_index];
  const std::vector<double>& energies = m_reader->GetEnergies();

  double previous_probability = 0.0;
  double m_rndenergy_probability = m_rnd.Rndm();
  for (size_t j = 0; j < energy_probability.size(); ++j) {
    if (energy_probability[j] >= m_rndenergy_probability) {
      // Apply linear interpolation to determine precise energy value
      double denominator = energy_probability[j] - previous_probability;
      double ratio = (denominator > s_epsilon) ? (energy_probability[j] - m_rndenergy_probability) / denominator : 0.0;

      double energy_low = energies[j];
      double energy_high = energies[j+1];
      m_energy_muon = energy_high - (energy_high - energy_low) * ratio;
      break;
    }
    previous_probability = energy_probability[j];
  }

  // 7. Final calculation of start/end coordinate
  double theta = TMath::DegToRad() * theta_deg;
  double tan_theta = TMath::Tan(theta);

  double dz_start = m_start_z - m_reference_z;
  m_start_x = m_reference_x + dz_start * tan_theta * cos_phi;
  m_start_y = m_reference_y + dz_start * tan_theta * sin_phi;
  double dz_end = m_end_z - m_reference_z;
  m_end_x = m_reference_x + dz_end * tan_theta * cos_phi;
  m_end_y = m_reference_y + dz_end * tan_theta * sin_phi;

  this->SetPoint(0, m_start_x, m_start_y, m_start_z);
  this->SetPoint(1, m_end_x, m_end_y, m_end_z);

  return true;
}

inline std::pair<double, double> Muon::CalculateEnergyLoss(const std::vector<Shape*>& shapes) {
  double energy_loss = 0.0;
  if (shapes.empty()) {
    if (g_debug_level>0) { Info("Muon::CalculateEnergyLoss", "No shape"); }
    return {m_energy_muon, energy_loss};
  }

  std::vector<TVector3> points;
  for (auto shape : shapes) {
    if (shape->IsHit(*this)) {
      std::vector<TVector3> hits = shape->HitPoints(*this);
      for (size_t i = 0; i < hits.size(); ++i) {
        if (g_debug_level>0) { Info("Muon::CalculateEnergyLoss", "hits [%zu]: (%.2f, %.2f, %.2f)", i, hits[i].X(), hits[i].Y(), hits[i].Z()); }
      }
      points.insert(points.end(), hits.begin(), hits.end());
    }
  }

  if (points.empty()) {
    if (g_debug_level>0) { Info("Muon::CalculateEnergyLoss", "No hit points"); }
    return {m_energy_muon, energy_loss};
  } else {
    if (g_debug_level>0) { Info("Muon::CalculateEnergyLoss", "Hit points size: %zu", points.size()); }
    // for (size_t i = 0; i < points.size(); ++i) {
    //   if (g_debug_level>0) { Info("Muon::CalculateEnergyLoss", "hit points [%zu]: (%.2f, %.2f, %.2f)", i, points[i].X(), points[i].Y(), points[i].Z()); }
    // }
  }

  // Sort points by Z descending (top to bottom)
  std::sort(points.begin(), points.end(), [](const TVector3& a, const TVector3& b) {
      return a.Z() > b.Z();
      });

  // Remove duplicates
  std::vector<TVector3> unique_points;
  unique_points.push_back(points[0]);
  for (size_t i = 1; i < points.size(); ++i) {
    if ((points[i] - unique_points.back()).Mag() > s_epsilon) {
      unique_points.push_back(points[i]);
    }
  }

  for (size_t i = 0; i < unique_points.size() - 1; ++i) {
    TVector3 segment_middle = (unique_points[i] + unique_points[i+1]) * 0.5;
    int highest_priority = Shape::s_invalid_priority;
    double selected_density = 0.0;
    bool found_material = false;

    if (g_debug_level>0) { Info("Muon::CalculateEnergyLoss", "highest priority: %d", highest_priority); }
    for (auto shape : shapes) {
      if (shape->IsInside(segment_middle)) {
        if (shape->GetPriority() < highest_priority) {
          highest_priority = shape->GetPriority();
          selected_density = shape->DensityAt(segment_middle);
          found_material = true;
        }
      }
    }

    if (g_debug_level>0) { 
      Info("Muon::CalculateEnergyLoss", "highest priority: %d", highest_priority);
      Info("Muon::CalculateEnergyLoss", "segment middle: (%.2f, %.2f, %.2f)", segment_middle.X(), segment_middle.Y(), segment_middle.Z());
      Info("Muon::CalculateEnergyLoss", "selected density: %.2f", selected_density);
    }
    if (found_material && selected_density > s_epsilon) {
      double current_muon_energy = m_energy_muon - energy_loss;
      if (current_muon_energy <= s_epsilon) break;
      double path_length = (unique_points[i] - unique_points[i+1]).Mag();
      energy_loss += CalculateSegmentEnergyLoss(selected_density, current_muon_energy, path_length);
    }
  }
  return {m_energy_muon, energy_loss};
}

inline double Muon::CalculateSegmentEnergyLoss(double density, double muon_energy, double path_length) {
  constexpr double kMuonMass = 105.658;   // [MeV]
  double current_energy = muon_energy;    // [MeV]
  double initial_energy = current_energy; // [MeV]

  double maximum_step = 1.0;  // [cm]
  double remaining_path = path_length;  // [cm]

  while (remaining_path > s_epsilon) {
    double step_length = TMath::Min(maximum_step, remaining_path);  // [cm]
                                                                    // Bethe-Bloch formula
    double dEdx = 1.88 + 0.077 * TMath::Log(current_energy / kMuonMass) + 3.9 * 1e-6 * current_energy;  // [MeV cm2/g]
    double loss = dEdx * density * step_length;  // [MeV cm2/g]*[g/cm3]*[cm] = [MeV]
    current_energy -= loss;
    remaining_path -= step_length;

    if (current_energy <= s_epsilon) {
      current_energy = 0.0;
      break;
    }
  }
  double segment_energy_loss = initial_energy - current_energy;
  return segment_energy_loss;
}

inline void Muon::AddHitMarker(const TVector3& point) {
  m_hit_markers->SetNextPoint(point.X(), point.Y(), point.Z());
}

inline void Muon::Draw(Option_t* option) {
  TPolyLine3D::Draw(option);

  TString opt = option;
  opt.ToLower();
  if (opt.Contains("m")) {
    if (m_hit_markers && m_hit_markers->GetN() > 0) m_hit_markers->Draw();
  }
}

inline void Muon::SetMarkerColor(Color_t color) {
  m_hit_markers->SetMarkerColor(color);
}

#endif

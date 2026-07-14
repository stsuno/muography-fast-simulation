#ifndef READ_EXPACS_H
#define READ_EXPACS_H

#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <sstream>

#include <TMath.h>
#include <TError.h>
#include <TSystem.h>

class ReadExpacs {
  public:
    ReadExpacs(double energy_threshold = 0.0, std::string mode_angle = "", std::string mode_energy = "left");

    // Accessors (Setters and Getters)
    void SetEnergyThreshold(double threshold);
    double GetEnergyThreshold() const { return m_energy_threshold; }
    // Basic Data
    const std::vector<double>& GetAngles() const { return m_angles; }
    const std::vector<double>& GetEnergies() const { return m_energies; }
    // Middle Points for Sampling (bin centers)
    const std::vector<double>& GetSamplingAngles() const { return m_sampling_angles; }
    const std::vector<double>& GetSamplingEnergies() const { return m_sampling_energies; }
    // Flux Matrices [angle][energy]
    const std::vector<std::vector<double>>& GetFluxP() const { return m_flux_p; }
    const std::vector<std::vector<double>>& GetFluxM() const { return m_flux_m; }
    const std::vector<std::vector<double>>& GetTresholdedFluxP() const { return m_thresholded_flux_p; }
    const std::vector<std::vector<double>>& GetThresholdedFluxM() const { return m_thresholded_flux_m; }
    const std::vector<std::vector<double>>& GetMiddleAngleFluxP() const { return m_middle_angle_flux_p; }
    const std::vector<std::vector<double>>& GetMiddleAngleFluxM() const { return m_middle_angle_flux_m; }
    // Integrated Flux Data
    const std::vector<double>& GetEnergyIntegratedFluxP() const { return m_energy_integrated_flux_p; }
    const std::vector<double>& GetEnergyIntegratedFluxM() const { return m_energy_integrated_flux_m; }
    const std::vector<double>& GetAngularFluxP() const { return m_angular_flux_p; }
    const std::vector<double>& GetAngularFluxM() const { return m_angular_flux_m; }
    double GetTotalFluxP() const { return m_total_flux_p; }
    double GetTotalFluxM() const { return m_total_flux_m; }
    // Cumulative Distribution Functions (CDFs)
    const std::vector<double>& GetAngleCDFP() const { return m_angle_cdf_p; }
    const std::vector<double>& GetAngleCDFM() const { return m_angle_cdf_m; }
    const std::vector<std::vector<double>>& GetEnergyCDFP() const { return m_energy_cdf_p; }
    const std::vector<std::vector<double>>& GetEnergyCDFM() const { return m_energy_cdf_m; }

  private:
    // Core Logic Functions
    bool LoadFile(const std::string& filename);
    void ProcessPhysicsData();
    void ApplyEnergyThreshold();
    void PrepareMiddleAngleFlux();
    void IntegrateFlux();
    void PrepareCDFs();

    // Utility Functions
    double GetRepresentativeValue(double v1, double v2, const std::string& mode) const;
    std::vector<double> Create1DCDF(const std::vector<double>& data, double total) const;
    std::vector<double> CreateEnergyCDF(size_t angle_index, double total, bool is_plus) const;

    // Member Variables
    double m_energy_threshold;
    std::string m_mode_angle;
    std::string m_mode_energy;

    std::vector<double> m_angles;
    std::vector<double> m_energies;
    std::vector<std::vector<double>> m_flux_p;
    std::vector<std::vector<double>> m_flux_m;

    std::vector<std::vector<double>> m_thresholded_flux_p;
    std::vector<std::vector<double>> m_thresholded_flux_m;
    std::vector<std::vector<double>> m_middle_angle_flux_p;
    std::vector<std::vector<double>> m_middle_angle_flux_m;

    std::vector<double> m_sampling_angles;
    std::vector<double> m_sampling_energies;
    std::vector<double> m_energy_steps;
    std::vector<double> m_energy_integrated_flux_p;
    std::vector<double> m_energy_integrated_flux_m;
    std::vector<double> m_angular_flux_p;
    std::vector<double> m_angular_flux_m;

    double m_total_flux_p = 0.0;
    double m_total_flux_m = 0.0;

    std::vector<double> m_angle_cdf_p;
    std::vector<double> m_angle_cdf_m;
    std::vector<std::vector<double>> m_energy_cdf_p;
    std::vector<std::vector<double>> m_energy_cdf_m;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor
inline ReadExpacs::ReadExpacs(double energy_threshold, std::string mode_angle, std::string mode_energy)
  : m_energy_threshold(energy_threshold), m_mode_angle(mode_angle), m_mode_energy(mode_energy) {
  if (LoadFile("data/input/expacs.output")) {
    ProcessPhysicsData();
  }
}

inline void ReadExpacs::SetEnergyThreshold(double threshold) {
  m_energy_threshold = threshold;
  Info("ReadExpacs::SetEnergyThreshold", "Recalculating for threshold: %.1f MeV", threshold);
  ProcessPhysicsData();
}

inline bool ReadExpacs::LoadFile(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    Error("ReadExpacs::LoadFile", "File not found: %s", filename.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }
  Info("ReadExpacs::LoadFile", "Read input file: %s", filename.c_str());

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::vector<double> values;
    double value;
    while (ss >> value) values.push_back(value);

    if (values.size() == 1) {
      m_angles.push_back(values[0]);
      m_flux_p.push_back({});
      m_flux_m.push_back({});
    } else if (values.size() >= 3) {
      if (m_angles.size() == 1) m_energies.push_back(values[0]);
      m_flux_p.back().push_back(values[1]);
      m_flux_m.back().push_back(values[2]);
    }
  }
  return true;
}

inline double ReadExpacs::GetRepresentativeValue(double v1, double v2, const std::string& mode) const {
  if (mode == "left") return v1;
  if (mode == "right") return v2;
  return (v1 + v2) * 0.5;
}

inline void ReadExpacs::ProcessPhysicsData() {
  ApplyEnergyThreshold();
  PrepareMiddleAngleFlux();
  IntegrateFlux();
  PrepareCDFs();
}

inline void ReadExpacs::ApplyEnergyThreshold() {
  Info("ReadExpacs::ApplyEnergyThreshold", "muon energy threshold: %.1f", m_energy_threshold);
  m_thresholded_flux_p.clear();
  m_thresholded_flux_m.clear();
  for (size_t i = 0; i < m_angles.size(); ++i) {
    std::vector<double> row_p, row_m;
    for (size_t j = 0; j < m_energies.size(); ++j) {
      if (m_energies[j] < m_energy_threshold) {
        row_p.push_back(0.0);
        row_m.push_back(0.0);
      } else {
        row_p.push_back(m_flux_p[i][j]);
        row_m.push_back(m_flux_m[i][j]);
      }
    }
    m_thresholded_flux_p.push_back(row_p);
    m_thresholded_flux_m.push_back(row_m);
  }
}

inline void ReadExpacs::PrepareMiddleAngleFlux() {
  m_middle_angle_flux_p.clear();
  m_middle_angle_flux_m.clear();
  for (size_t i = 0; i < m_angles.size() - 1; ++i) {
    std::vector<double> spec_p, spec_m;
    for (size_t j = 0; j < m_energies.size(); ++j) {
      spec_p.push_back(GetRepresentativeValue(m_thresholded_flux_p[i][j], m_thresholded_flux_p[i+1][j], m_mode_angle));
      spec_m.push_back(GetRepresentativeValue(m_thresholded_flux_m[i][j], m_thresholded_flux_m[i+1][j], m_mode_angle));
    }
    m_middle_angle_flux_p.push_back(spec_p);
    m_middle_angle_flux_m.push_back(spec_m);
  }
}

inline void ReadExpacs::IntegrateFlux() {
  m_sampling_angles.clear();
  m_sampling_energies.clear();
  m_energy_steps.clear();

  for (size_t i = 0; i < m_angles.size() - 1; ++i) {
    m_sampling_angles.push_back((m_angles[i] + m_angles[i+1]) * 0.5);
  }
  for (size_t j = 0; j < m_energies.size() - 1; ++j) {
    m_sampling_energies.push_back((m_energies[j] + m_energies[j+1]) * 0.5);
    m_energy_steps.push_back(m_energies[j+1] - m_energies[j]);
  }

  m_energy_integrated_flux_p.assign(m_sampling_angles.size(), 0.0);
  m_energy_integrated_flux_m.assign(m_sampling_angles.size(), 0.0);
  m_angular_flux_p.assign(m_sampling_angles.size(), 0.0);
  m_angular_flux_m.assign(m_sampling_angles.size(), 0.0);

  for (size_t i = 0; i < m_sampling_angles.size(); ++i) {
    for (size_t j = 0; j < m_energy_steps.size(); ++j) {
      double value_p = GetRepresentativeValue(m_middle_angle_flux_p[i][j], m_middle_angle_flux_p[i][j+1], m_mode_energy);
      double value_m = GetRepresentativeValue(m_middle_angle_flux_m[i][j], m_middle_angle_flux_m[i][j+1], m_mode_energy);
      m_energy_integrated_flux_p[i] += value_p * m_energy_steps[j];
      m_energy_integrated_flux_m[i] += value_m * m_energy_steps[j];
    }
    double theta_rad = TMath::DegToRad() * m_sampling_angles[i];
    double d_theta = TMath::DegToRad() * (m_angles[i+1] - m_angles[i]);
    // The EXPACS angular differential flux is a directional flux, i.e. per unit
    // area PERPENDICULAR to the particle direction, per steradian
    // [cm^-2 s^-1 sr^-1 MeV^-1] (T. Sato, PLOS ONE 11(8):e0160390, 2016).
    // Muons are generated on a HORIZONTAL plane with uniformly sampled (x, y),
    // so the crossing rate through that plane requires the projection factor
    // cos(theta) in addition to the solid-angle element 2*pi*sin(theta)*d_theta.
    double weight = 2.0 * TMath::Pi() * TMath::Sin(theta_rad) * TMath::Cos(theta_rad);
    m_angular_flux_p[i] = m_energy_integrated_flux_p[i] * weight * d_theta;
    m_angular_flux_m[i] = m_energy_integrated_flux_m[i] * weight * d_theta;
  }
  m_total_flux_p = std::accumulate(m_angular_flux_p.begin(), m_angular_flux_p.end(), 0.0);
  m_total_flux_m = std::accumulate(m_angular_flux_m.begin(), m_angular_flux_m.end(), 0.0);
  Info("ReadExpacs::IntegrateFlux", "Total Flux (plus) %f /cm2/s", m_total_flux_p);
  Info("ReadExpacs::IntegrateFlux", "Total Flux (minus) %f /cm2/s", m_total_flux_m);
}

inline void ReadExpacs::PrepareCDFs() {
  m_angle_cdf_p = Create1DCDF(m_angular_flux_p, m_total_flux_p);
  m_angle_cdf_m = Create1DCDF(m_angular_flux_m, m_total_flux_m);
  m_energy_cdf_p.clear(); m_energy_cdf_m.clear();
  for (size_t i = 0; i < m_sampling_angles.size(); ++i) {
    m_energy_cdf_p.push_back(CreateEnergyCDF(i, m_energy_integrated_flux_p[i], true));
    m_energy_cdf_m.push_back(CreateEnergyCDF(i, m_energy_integrated_flux_m[i], false));
  }
}

inline std::vector<double> ReadExpacs::Create1DCDF(const std::vector<double>& data, double total) const {
  std::vector<double> cdf;
  double current_sum = 0.0;
  for (double value : data) {
    if (total > 0) current_sum += value / total;
    cdf.push_back(current_sum);
  }
  return cdf;
}

inline std::vector<double> ReadExpacs::CreateEnergyCDF(size_t angle_index, double total, bool is_plus) const {
  std::vector<double> cdf;
  double current_sum = 0.0;
  const auto& target_spec = is_plus ? m_middle_angle_flux_p[angle_index] : m_middle_angle_flux_m[angle_index];
  for (size_t j = 0; j < m_energy_steps.size(); ++j) {
    double value = GetRepresentativeValue(target_spec[j], target_spec[j+1], m_mode_energy);
    if (total > 0) current_sum += (value * m_energy_steps[j]) / total;
    cdf.push_back(current_sum);
  }
  return cdf;
}

#endif

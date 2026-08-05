#ifndef SHAPE_H
#define SHAPE_H

#include <string>
#include <vector>

#include <TVector3.h>
#include <TObject.h>
#include <TAttLine.h>
#include <TAttFill.h>
#include <TError.h>
#include <TSystem.h>

// Forward declaration
class Muon;

class Shape : public TObject, public TAttLine, public TAttFill {
  public:
    // Sternheimer parameterization coefficients for the density-effect correction
    // delta(beta*gamma) used in the Bethe-Bloch formula (PDG Eq. 34.7).
    // Values for each material are provided by MaterialList in
    // Simulation/MaterialList.h, sourced from the PDG Atomic and Nuclear
    // Properties tables: https://pdg.lbl.gov/2026/AtomicNuclearProperties/
    struct SternheimerParameters {
      double m_a      = 0.0;  // coefficient of the intermediate-region term a*(x1-x)^k
      double m_k      = 0.0;  // exponent of the intermediate-region term (k = m_s)
      double m_x0     = 0.0;  // lower boundary of the intermediate region (x = log10(beta*gamma))
      double m_x1     = 0.0;  // upper boundary of the intermediate region
      double m_cbar   = 0.0;  // Cbar, the asymptotic offset: delta -> 2*ln(10)*x - Cbar
      double m_delta0 = 0.0;  // conductor correction below x0 (0 for nonconductors)
    };

    // Bundle of properties needed for muon transport (Bethe-Bloch), as
    // registered per material name in MaterialList (Simulation/MaterialList.h).
    // Materials for which only the density is known keep kUnset in the other
    // fields and a default-constructed (all-zero) m_sternheimer; SetMaterial()
    // below rejects those. Cbar > 0 for every real material, so m_cbar <= 0
    // identifies missing Sternheimer data (same convention as
    // Muon::CalculateDensityEffect).
    //
    // m_radiation_length is the only field not used by Bethe-Bloch: it feeds
    // the Highland multiple-Coulomb-scattering formula (PDG Eq. 34.16) in
    // Simulation/MuonTransport.h. It stays kUnset for materials whose X0 has
    // not been registered; SetMaterial() accepts those (scattering-free
    // transport still works), and MuonTransport reports them explicitly.
    struct Material {
      static constexpr double kUnset = -1.0;
      double m_density                = kUnset;  // [g/cm3]
      double m_z_over_a               = kUnset;  // <Z/A> [mol/g]
      double m_mean_excitation_energy = kUnset;  // I [MeV]
      SternheimerParameters m_sternheimer{};
      double m_radiation_length       = kUnset;  // X0 [g/cm2]
    };

    virtual ~Shape() {}

    // Mandatory functions for simulation
    // 1. Hit detection (returns true if the muon hits the shape)
    virtual bool IsHit(const Muon& muon) const = 0;
    // 2. Intersection points (returns a list of hit coordinates)
    virtual std::vector<TVector3> HitPoints(const Muon& muon) const = 0;
    // 3. Containment check (returns true if the point is inside the shape)
    virtual bool IsInside(const TVector3& point) const = 0;
    // 4. Density retrieval (returns density at a specific coordinate)
    virtual double DensityAt(const TVector3& point) const = 0;
    virtual double ZOverAAt(const TVector3& point) const = 0;
    virtual double MeanExcitationEnergyAt(const TVector3& point) const = 0;
    virtual SternheimerParameters SternheimerParametersAt(const TVector3& point) const = 0;

    // 5. Visualization
    virtual void Draw(Option_t* option = "") = 0;

    // 6. Geometry Operations
    virtual void MoveXYZ(double dx, double dy, double dz) = 0;
    virtual void RotatePhi(double phi) = 0;
    virtual void RotateTheta(double theta) = 0;
    virtual void MarkIntersections(Muon* muon) = 0;

    // Common properties
    virtual void SetDensity(double density) = 0;
    virtual double GetDensity() const = 0;
    virtual void SetPriority(int priority) = 0;
    virtual int GetPriority() const = 0;

    virtual void SetZOverA(double z_over_a) = 0;
    virtual double GetZOverA() const = 0;
    virtual void SetMeanExcitationEnergy(double mean_excitation_energy) = 0;
    virtual double GetMeanExcitationEnergy() const = 0;
    virtual void SetSternheimerParameters(const SternheimerParameters& parameters) = 0;
    virtual SternheimerParameters GetSternheimerParameters() const = 0;

    // Radiation length X0 [g/cm2], used by the multiple-scattering model in
    // Simulation/MuonTransport.h. Unlike the properties above these are not
    // pure virtual: storage lives in Shape, so existing subclasses keep
    // compiling unchanged. Subclasses whose material varies with position
    // (CompositeShape, Polyhedron) override RadiationLengthAt() to match the
    // behaviour of their DensityAt(). A value of 0 means "not set" and is
    // treated by MuonTransport as "do not scatter in this material".
    virtual void SetRadiationLength(double radiation_length) { m_radiation_length = radiation_length; }
    virtual double GetRadiationLength() const { return m_radiation_length; }
    virtual double RadiationLengthAt(const TVector3& /*point*/) const { return m_radiation_length; }

    // Convenience: sets density, <Z/A>, mean excitation energy, and Sternheimer
    // parameters at once from a Material bundle. Callers typically obtain one
    // via MaterialList::GetMaterial(name), e.g.:
    //   shape->SetMaterial(MaterialList::GetMaterial("ShieldingConcrete"));
    void SetMaterial(const Material& material);

    // Numerical Constants
    static constexpr double s_epsilon = 1e-9;  // Threshold to determine if a value is effectively zero
    static constexpr double s_tolerance = 1e-6;  // Tolerance for boundary and duplicate checks [cm] (10 nm)
    static constexpr int s_invalid_priority = 9999;  // Default value for invalid priority

  protected:
    double m_radiation_length = 0.0;  // X0 [g/cm2], 0 = not set
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline void Shape::SetMaterial(const Material& material) {
  if (material.m_z_over_a == Material::kUnset ||
      material.m_mean_excitation_energy == Material::kUnset ||
      material.m_sternheimer.m_cbar <= 0.0) {
    ::Error("Shape::SetMaterial", "Material lacks full Bethe-Bloch properties (density-only material?)");
    gSystem->Exit(EXIT_FAILURE);
  }
  SetDensity(material.m_density);
  SetZOverA(material.m_z_over_a);
  SetMeanExcitationEnergy(material.m_mean_excitation_energy);
  SetSternheimerParameters(material.m_sternheimer);
  // X0 is optional: an unregistered radiation length disables scattering in
  // this shape but must not block Bethe-Bloch transport, so no error here.
  SetRadiationLength(material.m_radiation_length == Material::kUnset ? 0.0
                                                                    : material.m_radiation_length);
}

#endif

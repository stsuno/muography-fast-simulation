#ifndef SHAPE_H
#define SHAPE_H

#include <vector>

#include <TVector3.h>
#include <TObject.h>
#include <TAttLine.h>
#include <TAttFill.h>

// Forward declaration
class Muon;

// Sternheimer parameterization coefficients for the density-effect correction
// delta(beta*gamma) used in the Bethe-Bloch formula (PDG Eq. 34.7).
// Values for each material are provided by GetSternheimerParameters() in
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

class Shape : public TObject, public TAttLine, public TAttFill {
  public:
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
    virtual double AtomicNumberAt(const TVector3& point) const = 0;
    virtual double AtomicMassAt(const TVector3& point) const = 0;
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

    virtual void SetAtomicNumber(double atomic_number) = 0;
    virtual double GetAtomicNumber() const = 0;
    virtual void SetAtomicMass(double atomic_mass) = 0;
    virtual double GetAtomicMass() const = 0;
    virtual void SetMeanExcitationEnergy(double mean_excitation_energy) = 0;
    virtual double GetMeanExcitationEnergy() const = 0;
    virtual void SetSternheimerParameters(const SternheimerParameters& parameters) = 0;
    virtual SternheimerParameters GetSternheimerParameters() const = 0;

    // Numerical Constants
    static constexpr double s_epsilon = 1e-9;  // Threshold to determine if a value is effectively zero
    static constexpr double s_tolerance = 1e-6;  // Tolerance for boundary and duplicate checks [cm] (10 nm)
    static constexpr int s_invalid_priority = 9999;  // Default value for invalid priority
};

#endif

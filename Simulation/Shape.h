#ifndef SHAPE_H
#define SHAPE_H

#include <vector>

#include <TVector3.h>

// Forward declaration
class Muon;

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

    // Numerical Constants
    static constexpr double s_epsilon = 1e-9;  // Threshold to determine if a value is effectively zero
    static constexpr double s_tolerance = 1e-6;  // Tolerance for boundary and duplicate checks [cm] (10 nm)
    static constexpr int s_invalid_priority = 9999;  // Default value for invalid priority
};

#endif

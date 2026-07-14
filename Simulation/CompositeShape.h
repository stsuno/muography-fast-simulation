#ifndef COMPOSITE_SHAPE_H
#define COMPOSITE_SHAPE_H

#include <vector>
#include <algorithm>

#include <TVector3.h>
#include <TMath.h>

#include "Simulation/Shape.h"

class CompositeShape : public Shape {
  public:
    CompositeShape() = default;
    virtual ~CompositeShape();

    // Manage child shapes
    void AddShape(Shape* shape);

    // Shape Interface Implementation
    virtual bool IsHit(const Muon& muon) const override;
    virtual std::vector<TVector3> HitPoints(const Muon& muon) const override;
    virtual bool IsInside(const TVector3& point) const override;
    virtual double DensityAt(const TVector3& point) const override;
    virtual double AtomicNumberAt(const TVector3& point) const override;
    virtual double AtomicMassAt(const TVector3& point) const override;
    virtual double MeanExcitationEnergyAt(const TVector3& point) const override;
    virtual SternheimerParameters SternheimerParametersAt(const TVector3& point) const override;


    // Visualization
    virtual void Draw(Option_t* option = "") override;
    virtual void SetLineColor(Color_t color) override;
    virtual void SetLineWidth(Width_t width) override;
    virtual void SetLineStyle(Style_t style) override;
    virtual void SetFillColor(Color_t color) override;
    virtual void SetFillStyle(Style_t style) override;

    // Geometry Operations
    void MoveXYZ(double dx, double dy, double dz) override;
    void RotatePhi(double phi) override;
    void RotateTheta(double theta) override;
    void MarkIntersections(Muon* muon) override;

    // Attributes (Using the first shape as a representative)
    virtual void SetDensity(double density) override;
    virtual double GetDensity() const override { return m_density; }
    virtual void SetPriority(int priority) override;
    virtual int GetPriority() const override { return m_priority; }

    virtual void SetAtomicNumber(double atomic_number) override;
    virtual double GetAtomicNumber() const override { return m_atomic_number; }
    virtual void SetAtomicMass(double atomic_mass) override;
    virtual double GetAtomicMass() const override { return m_atomic_mass; }
    virtual void SetMeanExcitationEnergy(double mean_excitation_energy) override;
    virtual double GetMeanExcitationEnergy() const override { return m_mean_excitation_energy; }
    virtual void SetSternheimerParameters(const SternheimerParameters& parameters) override;
    virtual SternheimerParameters GetSternheimerParameters() const override { return m_sternheimer_parameters; }

  private:
    // List of shapes, ordered by priority (index 0 is highest)
    std::vector<Shape*> m_shapes;

    double m_density = 0.0;
    int m_priority = s_invalid_priority;

    double m_atomic_number = 0.0;
    double m_atomic_mass = 0.0;
    double m_mean_excitation_energy = 0.0;
    SternheimerParameters m_sternheimer_parameters;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline CompositeShape::~CompositeShape() {
  for (auto shape : m_shapes) {
    if (shape) delete shape;
  }
  m_shapes.clear();
}

inline void CompositeShape::AddShape(Shape* shape) {
  if (shape) m_shapes.push_back(shape);
}

inline bool CompositeShape::IsHit(const Muon& muon) const {
  for (auto shape : m_shapes) {
    if (shape->IsHit(muon)) return true;
  }
  return false;
}

inline std::vector<TVector3> CompositeShape::HitPoints(const Muon& muon) const {
  std::vector<TVector3> all_hit_points;
  for (auto shape : m_shapes) {
    std::vector<TVector3> hit_points = shape->HitPoints(muon);
    all_hit_points.insert(all_hit_points.end(), hit_points.begin(),  hit_points.end());
  }

  std::vector<TVector3> unique_hit_points;
  for (auto& point : all_hit_points) {
    auto it = std::find_if(unique_hit_points.begin(), unique_hit_points.end(),
        [&](const TVector3& hit) {return (hit - point).Mag() < s_tolerance; });
    if (it == unique_hit_points.end()) {
      unique_hit_points.push_back(point);
    }
  }

  // Sort points by Z coordinate (Descending: from Top to Bottom)
  std::sort(unique_hit_points.begin(), unique_hit_points.end(),
      [](const TVector3& a, const TVector3& b) { return a.Z() > b.Z(); });

  return unique_hit_points;
}

inline bool CompositeShape::IsInside(const TVector3& point) const {
  for (auto shape : m_shapes) {
    if (shape->IsInside(point)) return true;
  }
  return false;
}

inline double CompositeShape::DensityAt(const TVector3& point) const {
  // Returns the density of the first shape that contains the point
  for (auto shape : m_shapes) {
    if (shape->IsInside(point)) {
      return shape->GetDensity();
    }
  }
  return 0.0;
}

inline double CompositeShape::AtomicNumberAt(const TVector3& point) const {
  // Returns the density of the first shape that contains the point
  for (auto shape : m_shapes) {
    if (shape->IsInside(point)) {
      return shape->GetAtomicNumber();
    }
  }
  return 0.0;
}

inline double CompositeShape::AtomicMassAt(const TVector3& point) const {
  // Returns the density of the first shape that contains the point
  for (auto shape : m_shapes) {
    if (shape->IsInside(point)) {
      return shape->GetAtomicMass();
    }
  }
  return 0.0;
}

inline double CompositeShape::MeanExcitationEnergyAt(const TVector3& point) const {
  // Returns the density of the first shape that contains the point
  for (auto shape : m_shapes) {
    if (shape->IsInside(point)) {
      return shape->GetMeanExcitationEnergy();
    }
  }
  return 0.0;
}

inline SternheimerParameters CompositeShape::SternheimerParametersAt(const TVector3& point) const {
  // Returns the parameters of the first shape that contains the point
  for (auto shape : m_shapes) {
    if (shape->IsInside(point)) {
      return shape->GetSternheimerParameters();
    }
  }
  return SternheimerParameters{};
}

inline void CompositeShape::MoveXYZ(double dx, double dy, double dz) {
  for (auto shape : m_shapes) shape->MoveXYZ(dx, dy, dz);  
}

inline void CompositeShape::RotatePhi(double phi) {
  for (auto shape : m_shapes) shape->RotatePhi(phi);
}

inline void CompositeShape::RotateTheta(double theta) {
  for (auto shape : m_shapes) shape->RotateTheta(theta);
}

inline void CompositeShape::MarkIntersections(Muon* muon) {
  for (auto shape : m_shapes) shape->MarkIntersections(muon);
}

inline void CompositeShape::Draw(Option_t* option) {
  for (auto shape : m_shapes) shape->Draw(option);
}

inline void CompositeShape::SetDensity(double density) {
  m_density = density;
  for (auto shape : m_shapes) shape->SetDensity(density);
}

inline void CompositeShape::SetPriority(int priority) {
  m_priority = priority;
  for (auto shape : m_shapes) shape->SetPriority(priority);
}

inline void CompositeShape::SetAtomicNumber(double atomic_number) {
  m_atomic_number = atomic_number;
  for (auto shape : m_shapes) shape->SetAtomicNumber(atomic_number);
}

inline void CompositeShape::SetAtomicMass(double atomic_mass) {
  m_atomic_mass = atomic_mass;
  for (auto shape : m_shapes) shape->SetAtomicMass(atomic_mass);
}

inline void CompositeShape::SetMeanExcitationEnergy(double mean_excitation_energy) {
  m_mean_excitation_energy = mean_excitation_energy;
  for (auto shape : m_shapes) shape->SetMeanExcitationEnergy(mean_excitation_energy);
}

inline void CompositeShape::SetSternheimerParameters(const SternheimerParameters& parameters) {
  m_sternheimer_parameters = parameters;
  for (auto shape : m_shapes) shape->SetSternheimerParameters(parameters);
}

inline void CompositeShape::SetLineColor(Color_t color) {
  TAttLine::SetLineColor(color);
  for (auto shape : m_shapes) {
    if (shape) shape->SetLineColor(color);
  }
}

inline void CompositeShape::SetLineWidth(Width_t width) {
  TAttLine::SetLineWidth(width);
  for (auto shape : m_shapes) {
    if (shape) shape->SetLineWidth(width);
  }
}

inline void CompositeShape::SetLineStyle(Style_t style) {
  TAttLine::SetLineStyle(style);
  for (auto shape : m_shapes) {
    if (shape) shape->SetLineStyle(style);
  }
}

inline void CompositeShape::SetFillColor(Color_t color) {
  TAttFill::SetFillColor(color);
  for (auto shape : m_shapes) {
    if (shape) shape->SetFillColor(color);
  }
}

inline void CompositeShape::SetFillStyle(Style_t style) {
  TAttFill::SetFillStyle(style);
  for (auto shape : m_shapes) {
    if (shape) shape->SetFillStyle(style);
  }
}

#endif

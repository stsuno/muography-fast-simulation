#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <vector>
#include <optional>
#include <algorithm>

#include <TPolyLine3D.h>
#include <TVector3.h>
#include <TMath.h>

#include "Simulation/Muon.h"
#include "Simulation/Shape.h"

// Structure to hold coefficients of the plane equation: ax + by + cz + d = 0
struct SurfaceCoefficient {
  double m_a;  // coefficient for x
  double m_b;  // coefficient for y
  double m_c;  // coefficient for z
  double m_d;  // constant term
};

class Polyhedron : public Shape {
  public:
    Polyhedron() = default;
    virtual ~Polyhedron();

    // Shape Interface Implemetntation
    bool IsHit(const Muon& muon) const override;
    std::vector<TVector3> HitPoints(const Muon& muon) const override;
    bool IsInside(const TVector3& point) const override;
    double DensityAt(const TVector3& point) const override;
    double AtomicNumberAt(const TVector3& point) const override;
    double AtomicMassAt(const TVector3& point) const override;
    double MeanExcitationEnergyAt(const TVector3& point) const override;
    SternheimerParameters SternheimerParametersAt(const TVector3& point) const override;


    // Visualization
    void Draw(Option_t* option = "") override;
    void Fill(Option_t* option = "");
    void MarkIntersections(Muon* muon) override;

    // Geometry Operations
    void MoveXYZ(double dx, double dy, double dz) override;
    void RotatePhi(double phi) override;
    void RotateTheta(double theta) override;

    // Accessors (Setters and Getters)
    void SetDensity(double density) override { m_density = density; }
    double GetDensity() const override { return m_density; }
    void SetPriority(int priority) override { m_priority = priority; }
    int GetPriority() const override { return m_priority; }

    double GetAtomicNumber() const override { return m_atomic_number; }
    void SetAtomicNumber(double atomic_number) override { m_atomic_number = atomic_number; }
    double GetAtomicMass() const override { return m_atomic_mass; }
    void SetAtomicMass(double atomic_mass) override { m_atomic_mass = atomic_mass; }
    double GetMeanExcitationEnergy() const override { return m_mean_excitation_energy; }
    void SetMeanExcitationEnergy(double mean_excitation_energy) override { m_mean_excitation_energy = mean_excitation_energy; }
    SternheimerParameters GetSternheimerParameters() const override { return m_sternheimer_parameters; }
    void SetSternheimerParameters(const SternheimerParameters& parameters) override { m_sternheimer_parameters = parameters; }

    TVector3 GetCenter() const { return m_center; }

    std::vector<TVector3> GetVertices() const { return m_vertices; }
    void UpdateVertex(int index, TVector3 vertex) { m_vertices.at(index)=vertex; }

    void ComputeCenter();
    void ComputeSurfaceEquation();

  protected:
    // Geometric Calculations
    std::optional<TVector3> GetIntersection(const Muon& muon, int face_index) const;

    // Member variables
    std::vector<TVector3> m_vertices;
    std::vector<std::vector<int>> m_faces;
    std::vector<bool> m_faces_flag;
    TVector3 m_center;
    std::vector<SurfaceCoefficient> m_coefficients;

    std::vector<TPolyLine3D*> m_wireframes {};
    std::vector<TPolyLine3D*> m_fill_lines {};

    double m_density = 0.0;
    int m_priority = s_invalid_priority;

    double m_atomic_number = 0.0;
    double m_atomic_mass = 0.0;
    double m_mean_excitation_energy = 0.0;
    SternheimerParameters m_sternheimer_parameters;

    void SetFace(const std::vector<std::vector<int>>& faces);
    void SetWireframe(int face_index);
    void SetFillLine(int face_index);
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline Polyhedron::~Polyhedron() {
  for (auto line : m_wireframes) {
    if (line) delete line;
  }
  for (auto line : m_fill_lines) {
    if (line) delete line;
  }
  m_wireframes.clear();
  m_fill_lines.clear();
}

inline void Polyhedron::SetFace(const std::vector<std::vector<int>>& faces) {
  m_faces = faces;
  m_faces_flag.assign(m_faces.size(), true);
}

inline void Polyhedron::ComputeCenter() {
  if (m_vertices.empty()) return;
  TVector3 sum(0, 0, 0);
  for (const auto& vertex : m_vertices) sum += vertex;
  m_center = sum * (1.0 / static_cast<double>(m_vertices.size()));
}

inline void Polyhedron::ComputeSurfaceEquation() {
  m_coefficients.clear();
  for (const auto& face : m_faces) {
    if (face.size() < 3) continue;
    // Get the first three vertices of the face
    TVector3 vertex0 = m_vertices[face[0]];
    TVector3 vertex1 = m_vertices[face[1]];
    TVector3 vertex2 = m_vertices[face[2]];
    // Define two vectors on the plane
    TVector3 vector1 = vertex1 - vertex0;
    TVector3 vector2 = vertex2 - vertex0;
    // Calculate the normal vector using cross product
    TVector3 normal = vector1.Cross(vector2).Unit();  // Normalize to [a, b, c]
    double d = -normal.Dot(vertex0);

    m_coefficients.push_back({normal.X(), normal.Y(), normal.Z(), d});
    // Info("Polyhedron::ComputeSurfaceEquation", "coefficient: (%.2f, %.2f, %.2f, %.2f)", normal.X(), normal.Y(), normal.Z(), d);
  }
}

inline void Polyhedron::SetWireframe(int face_index) {
  const std::vector<int>& face = m_faces[face_index];
  int n_points = face.size();
  TPolyLine3D* line = new TPolyLine3D(n_points + 1);

  for (int i = 0; i < n_points; ++i) {
    const TVector3& vertex = m_vertices[face[i]];
    line->SetPoint(i, vertex.X(), vertex.Y(), vertex.Z());
  }
  const TVector3& vertex0 = m_vertices[face[0]];
  line->SetPoint(n_points, vertex0.X(), vertex0.Y(), vertex0.Z());

  m_wireframes.push_back(line);
}

inline void Polyhedron::SetFillLine(int face_index) {
  const std::vector<int>& face = m_faces[face_index];
  int n_vertices = face.size();

  if (n_vertices == 3) {
    // Triangle: Interpolate between two edges (v0-v2 and v1-v2).
    const TVector3& vertex0 = m_vertices[face[0]];
    const TVector3& vertex1 = m_vertices[face[1]];
    const TVector3& vertex2 = m_vertices[face[2]];

    int grid = static_cast<int>((vertex0 - vertex2).Mag());
    if (grid < 3) grid = 3;

    TPolyLine3D* line = new TPolyLine3D((grid + 1) * 2);
    for (int i = 0; i <= grid; ++i) {
      double ratio = static_cast<double>(2 * i + 1) / ((grid + 1) * 2);

      // Create lines that are roughly parallel to the vase v0-v1
      TVector3 point0 = (1.0 - ratio) * vertex0 + ratio * vertex2;
      TVector3 point1 = (1.0 - ratio) * vertex1 + ratio * vertex2;

      line->SetPoint(2*i  , point0.X(), point0.Y(), point0.Z());
      line->SetPoint(2*i+1, point1.X(), point1.Y(), point1.Z());
    }
    m_fill_lines.push_back(line);

  } else if (n_vertices == 4) {
    // Quadrilateral: Interpolate between two opposite edges (v0-v1 and v3-v2).
    const TVector3& vertex0 = m_vertices[face[0]];
    const TVector3& vertex1 = m_vertices[face[1]];
    const TVector3& vertex2 = m_vertices[face[2]];
    const TVector3& vertex3 = m_vertices[face[3]];

    int grid = static_cast<int>((vertex0 - vertex1).Mag());
    if (grid < 3) grid = 3;

    TPolyLine3D* line = new TPolyLine3D((grid + 1) * 2);
    for (int i = 0; i <= grid; ++i) {
      double ratio = static_cast<double>(2 * i + 1) / ((grid + 1) * 2);

      // Create lines that are roughly parallel to the vase v1-v2
      TVector3 point0 = (1.0 - ratio) * vertex0 + ratio * vertex1;
      TVector3 point1 = (1.0 - ratio) * vertex3 + ratio * vertex2;

      line->SetPoint(2*i  , point0.X(), point0.Y(), point0.Z());
      line->SetPoint(2*i+1, point1.X(), point1.Y(), point1.Z());
    }
    m_fill_lines.push_back(line);
  }
}

inline void Polyhedron::Draw(Option_t* option) {
  // Clear previous wireframes
  for (auto line : m_wireframes) {
    if (line) delete line;
  }
  m_wireframes.clear();

  TString opt = option;
  opt.ToLower();
  if (opt.Contains("f")) this->Fill(option);

  for (size_t i = 0; i < m_faces.size(); ++i) {
    if (!m_faces_flag[i]) continue;
    SetWireframe(static_cast<int>(i));
  }

  for (auto line : m_wireframes) {
    if (!line) return;
    line->SetLineColor(this->GetLineColor());
    line->SetLineWidth(this->GetLineWidth());
    line->SetLineStyle(this->GetLineStyle());
    line->Draw(option);
  }
}

inline void Polyhedron::Fill(Option_t* option) {
  // Clear previous fill lines
  for (auto line : m_fill_lines) delete line;
  m_fill_lines.clear();

  TString opt = option;
  opt.ToLower();

  for (size_t i = 0; i < m_faces.size(); ++i) {
    if (!m_faces_flag[i]) continue;
    SetFillLine(static_cast<int>(i));
  }

  for (auto line : m_fill_lines) {
    if (!line) return;
    // Use FillColor as the line color for the mesh faces
    line->SetLineColorAlpha(this->GetFillColor(), 0.02);  // Clear
    if (opt.Contains("d")) line->SetLineColorAlpha(this->GetFillColor(), 1.0);  // Dense
                                                                                // Redirect FillStyle to LineStyle to control mesh appearance
    line->SetLineStyle(this->GetFillStyle());
    line->SetLineWidth(1);
    line->Draw(option);
  }
}

inline void Polyhedron::MarkIntersections(Muon* muon) {
  for (size_t i = 0; i < m_faces.size(); ++i) {
    std::optional<TVector3> intersection = GetIntersection(*muon, i);
    if (intersection.has_value()) muon->AddHitMarker(intersection.value());
  }
}

// Intersection calculation between a muon track and a specific face (OBB)
// Returns the intersection point if found, otherwise std::nullopt
inline std::optional<TVector3> Polyhedron::GetIntersection(const Muon& muon, int face_index) const {
  TVector3 start_point = muon.GetStartPoint();
  TVector3 direction = muon.GetEndPoint() - start_point;

  const SurfaceCoefficient& coefficient = m_coefficients[face_index];
  TVector3 normal(coefficient.m_a, coefficient.m_b, coefficient.m_c);
  double d = coefficient.m_d;

  double denominator = normal.Dot(direction);
  if (TMath::Abs(denominator) < s_epsilon) return std::nullopt;

  double t = -(normal.Dot(start_point) + d) / denominator;
  TVector3 intersection = start_point + t * direction;
  // Info("Polyhedron::GetIntersection", "Intersection: (%.2f, %.2f, %.2f)", intersection.X(), intersection.Y(), intersection.Z());

  // Check if intersection is within face boundary
  const std::vector<int>& face_indices = m_faces[face_index];
  bool is_all_positive = true;
  bool is_all_negative = true;

  for (size_t i = 0; i < face_indices.size(); ++i) {
    const TVector3& vertex_current = m_vertices[face_indices[i]];
    const TVector3& vertex_next = m_vertices[face_indices[(i + 1) % face_indices.size()]];
    TVector3 vector_edge = vertex_next - vertex_current;
    TVector3 vector_to_intersection = intersection - vertex_current;
    TVector3 cross_product = vector_edge.Cross(vector_to_intersection);
    double dot_product = normal.Dot(cross_product);
    // Precision-aware boundary check using kTolerance
    if (dot_product < -s_tolerance) is_all_positive = false;
    if (dot_product > s_tolerance) is_all_negative = false;
  }

  // If the point is consistently on one side of all edges, it's inside the face boundary
  if (is_all_positive || is_all_negative) return intersection;

  return std::nullopt;
}

// Check if the muon hits any of the polyhedron's faces
inline bool Polyhedron::IsHit(const Muon& muon) const {
  for (size_t i = 0; i < m_faces.size(); ++i) {
    if (!m_faces_flag[i]) continue;
    if (GetIntersection(muon, i).has_value()) return true;
  }
  return false;
}

// Get the entry and exit points of the muon through the polyhedron
inline std::vector<TVector3> Polyhedron::HitPoints(const Muon& muon) const {
  std::vector<TVector3> hit_points;
  for (size_t i = 0; i < m_faces.size(); ++i) {
    if (!m_faces_flag[i]) continue;

    // Calculate intersection for each face
    std::optional<TVector3> intersection = GetIntersection(muon, i);
    if (intersection.has_value()) {
      hit_points.push_back(intersection.value());
      // Info("Polyhedron::HitPoints", "hit point: (%.2f, %.2f, %.2f)", intersection.value().X(), intersection.value().Y(), intersection.value().Z());
    }
  }
  if (hit_points.empty()) return {};
  // Info("Polyhedron::Hitpoints", "hit_points size: %zu", hit_points.size());

  // for (size_t i = 0; i < hit_points.size(); ++i)  Info("Polyhedron::Hitpoints", "hit_points [%zu]: (%.2f, %.2f, %.2f)", i, hit_points[i].X(), hit_points[i].Y(), hit_points[i].Z());

  // Remove duplicates using a distance-based tolerance
  std::vector<TVector3> unique_hit_points;
  for (auto& point : hit_points) {
    auto it = std::find_if(unique_hit_points.begin(), unique_hit_points.end(),
        [&](const TVector3& hit) {return (hit - point).Mag() < s_tolerance; });
    if (it == unique_hit_points.end()) {
      unique_hit_points.push_back(point);
    }
  }
  // Info("Polyhedron::Hitpoints", "unique_hit_points size: %zu", unique_hit_points.size());

  // Sort points by Z coordinate (Descending: from Top to Bottom)
  std::sort(unique_hit_points.begin(), unique_hit_points.end(),
      [](const TVector3& a, const TVector3& b) { return a.Z() > b.Z(); });

  // for (size_t i = 0; i < unique_hit_points.size(); ++i) Info("Polyhedron::Hitpoints", "unique_hit_points [%zu]: (%.2f, %.2f, %.2f)", i, unique_hit_points[i].X(), unique_hit_points[i].Y(), unique_hit_points[i].Z());
  return unique_hit_points;
}

// Check if a point is inside the polyhedron
inline bool Polyhedron::IsInside(const TVector3& point) const {
  if (m_coefficients.empty()) return false;
  for (const auto& coefficient : m_coefficients) {
    // Calculate f = ax + by + cz + d for the given point
    double point_value = coefficient.m_a * point.X() + coefficient.m_b * point.Y() +
      coefficient.m_c * point.Z() + coefficient.m_d;
    // Calculate f = ax + by + cz + d for the center of the polyhedron
    double center_value = coefficient.m_a * m_center.X() + coefficient.m_b * m_center.Y() +
      coefficient.m_c * m_center.Z() + coefficient.m_d;
    // If point_value * center_value < 0, the point and center are on opposite sides of the plane.
    // Include points on the boundary by applying s_epsilon
    // Info("Polyhedron::IsInside", "point * center: %.2f", point_value * center_value);
    if (point_value * center_value < -s_epsilon) return false;
  }
  //Point is on the same side as the center for all faces
  return true;
}

inline double Polyhedron::DensityAt(const TVector3& point) const {
  return IsInside(point) ? m_density : 0.0;
}

inline double Polyhedron::AtomicNumberAt(const TVector3& point) const {
  return IsInside(point) ? m_atomic_number : 0.0;
}

inline double Polyhedron::AtomicMassAt(const TVector3& point) const {
  return IsInside(point) ? m_atomic_mass : 0.0;
}

inline double Polyhedron::MeanExcitationEnergyAt(const TVector3& point) const {
  return IsInside(point) ? m_mean_excitation_energy : 0.0;
}

inline SternheimerParameters Polyhedron::SternheimerParametersAt(const TVector3& point) const {
  return IsInside(point) ? m_sternheimer_parameters : SternheimerParameters{};
}

inline void Polyhedron::MoveXYZ(double dx, double dy, double dz) {
  TVector3 displacement(dx, dy, dz);
  for (auto& vertex : m_vertices) {
    vertex += displacement;
  }
  ComputeCenter();
  ComputeSurfaceEquation();
}

inline void Polyhedron::RotatePhi(double phi) {
  for (auto& vertex : m_vertices) {
    vertex -= m_center;
    vertex.RotateZ(phi);
    vertex += m_center;
  }
  ComputeCenter();
  ComputeSurfaceEquation();
}

inline void Polyhedron::RotateTheta(double theta) {
  for (auto& vertex : m_vertices) {
    vertex -= m_center;
    vertex.RotateX(theta);
    vertex += m_center;
  }
  ComputeCenter();
  ComputeSurfaceEquation();
}

#endif

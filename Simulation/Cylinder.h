#ifndef CYLINDER_H
#define CYLINDER_H

#include <TMath.h>

#include "Simulation/Polyhedron.h"

class Cylinder : public Polyhedron {
  public:
    Cylinder(double cx, double cy, double cz, double length, double radius);
    virtual ~Cylinder() = default;

  private:
    void SetVertex(double length, double radius);
    int m_num_edges = 12;  // n-gon
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline Cylinder::Cylinder(double cx, double cy, double cz, double length, double radius){

  std::vector<std::vector<int>> faces;
  // Bottom face (lower cap)
  std::vector<int> bottom;
  for (int i = 0; i < m_num_edges; ++i) bottom.push_back(i);
  faces.push_back(bottom);
  // Top face (upper cap)
  std::vector<int> top;
  for (int i = 0; i < m_num_edges; ++i) top.push_back(m_num_edges + i);
  faces.push_back(top);
  // Side faces (quadrilaterals connecting bottom and top edges)
  for (int i = 0; i < m_num_edges; ++i) {
    int next = (i + 1) % m_num_edges;
    // Indicies: bottom[i], bottom[next], top[next], top[i]
    faces.push_back({i, next, m_num_edges + next, m_num_edges + i});
  }

  SetFace(faces);
  m_faces_flag[0] = true;  // bottom
  m_faces_flag[1] = true;  // top

  m_vertices.resize(m_num_edges * 2);
  SetVertex(length, radius);
  ComputeCenter();
  MoveXYZ(cx - m_center.X(), cy - m_center.Y(), cz - m_center.Z());
}

inline void Cylinder::SetVertex(double length, double radius) {
  const double d_phi = TMath::TwoPi() / m_num_edges;
  for (int i = 0; i < m_num_edges; ++i) {
    double phi = i * d_phi;
    double x = radius * TMath::Cos(phi);
    double y = radius * TMath::Sin(phi);

    // Bottom vertices: index 0 to n_gon-1
    m_vertices[i].SetXYZ(x, y, 0.0);
    // Top vertices: index n_gon to 2*n_gon-1
    m_vertices[m_num_edges + i].SetXYZ(x, y, length);
  }
}

#endif

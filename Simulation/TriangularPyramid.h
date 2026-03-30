#ifndef TRIANGULAR_PYRAMID_H
#define TRIANGULAR_PYRAMID_H

#include "Simulation/Polyhedron.h"


class TriangularPyramid : public Polyhedron {
public:
  TriangularPyramid(double cx, double cy, double cz,
		    double width, double depth, double height,
		    bool is_right_angle = false);
  virtual ~TriangularPyramid() = default;

private:
  void SetVertex(double width, double height, double length, bool is_right_angle);
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline TriangularPyramid::TriangularPyramid(double cx, double cy, double cz,
					    double width, double depth, double height,
					    bool is_right_angle) {

  SetFace({       //                            *  right-angled / isoscales *
      {0, 1, 2},  // bottom      3              *                           *
      {0, 1, 3},  // front              z  y    *   y       2         2     *
      {1, 2, 3},  // right       2      | /     *   |       |\       / \    * 
      {2, 0, 3}   // left       / \     |/      *   |       | \     /   \   *
    });           //           0---1    o---x   *   o--x (3)0--1   0-(3)-1  *

  m_vertices.resize(4);
  SetVertex(width, depth, height, is_right_angle);
  ComputeCenter();
  MoveXYZ(cx - m_center.X(), cy - m_center.Y(), cz - m_center.Z());
}

inline void TriangularPyramid::SetVertex(double width, double depth, double height, bool is_right_angle) {
  if (is_right_angle) {
    // Right-angled triangle
    m_vertices[0].SetXYZ(  0.0,   0.0,    0.0);
    m_vertices[1].SetXYZ(width,   0.0,    0.0);
    m_vertices[2].SetXYZ(  0.0, depth,    0.0);
    m_vertices[3].SetXYZ(  0.0,   0.0, height);
  } else {
    // Isoscales triangle
    m_vertices[0].SetXYZ(        0.0,   0.0,    0.0);
    m_vertices[1].SetXYZ(      width,   0.0,    0.0);
    m_vertices[2].SetXYZ(0.5 * width, depth,    0.0);
    m_vertices[3].SetXYZ(0.5 * width,   0.0, height);
  }
}

#endif

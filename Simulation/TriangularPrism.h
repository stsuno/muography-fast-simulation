#ifndef TRIANGULAR_PRISM_H
#define TRIANGULAR_PRISM_H

#include "Simulation/Polyhedron.h"


class TriangularPrism : public Polyhedron {
public:
  TriangularPrism(double cx, double cy, double cz,
		  double width, double depth, double height,
		  bool is_right_angle = false);
  virtual ~TriangularPrism() = default;

private:
  void SetVertex(double width, double height, double length, bool is_right_angle);
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline TriangularPrism::TriangularPrism(double cx, double cy, double cz,
					double width, double depth, double height,
					bool is_right_angle) {

  SetFace({          //             5    x: width  *                          *
      {0, 1, 2},     // bottom     /|\   y: depth  * right-angled / isoscales *
      {3, 4, 4},     // top       3---4  z: height *                          *
      {0, 1, 4, 3},  // front     | | |    z  y    *   y      2         2     *
      {1, 2, 5, 4},  // right     |(2)|    | /     *   |      |\       / \    *
      {2, 0, 3, 5},  // left      |/ \|    |/      *   |      | \     /   \   *
    });              // right     0---1    o---x   *   o--x   0--1   0-----1  *

  m_vertices.resize(6);
  SetVertex(width, depth, height, is_right_angle);
  ComputeCenter();
  MoveXYZ(cx - m_center.X(), cy - m_center.Y(), cz - m_center.Z());
}

inline void TriangularPrism::SetVertex(double width, double depth, double height, bool is_right_angle) {
  if (is_right_angle) {
    // Right-angled triangle
    m_vertices[0].SetXYZ(  0.0,   0.0,    0.0);
    m_vertices[1].SetXYZ(width,   0.0,    0.0);
    m_vertices[2].SetXYZ(  0.0, depth,    0.0);
    m_vertices[3].SetXYZ(  0.0,   0.0, height);
    m_vertices[4].SetXYZ(width,   0.0, height);
    m_vertices[5].SetXYZ(  0.0, depth, height);
  } else {
    // Isoscales triangle
    m_vertices[0].SetXYZ(        0.0,   0.0,    0.0);
    m_vertices[1].SetXYZ(      width,   0.0,    0.0);
    m_vertices[2].SetXYZ(0.5 * width, depth,    0.0);
    m_vertices[3].SetXYZ(        0.0,   0.0, height);
    m_vertices[4].SetXYZ(      width,   0.0, height);
    m_vertices[5].SetXYZ(0.5 * width, depth, height);
  }
}

#endif

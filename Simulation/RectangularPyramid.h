#ifndef RECTANGULAR_PYRAMID_H
#define RECTANGULAR_PYRAMID_H

#include "Simulation/Polyhedron.h"


class RectangularPyramid : public Polyhedron {
  public:
    RectangularPyramid(double cx, double cy, double cz, double wx, double wy, double wz);
    virtual ~RectangularPyramid() = default;

  private:
    void SetVertex(double wx, double wy, double wz);
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline RectangularPyramid::RectangularPyramid(double cx, double cy, double cz, double wx, double wy, double wz) {

  SetFace({
      {0, 1, 2, 3},  // bottom         4
      {0, 1, 4},     // front                  z  y
      {1, 2, 4},     // right       3------2   | /
      {2, 3, 4},     // back       /      /    |/
      {3, 0, 4}      // back      0------1     o---x
    });

  m_vertices.resize(5);
  SetVertex(wx, wy, wz);
  ComputeCenter();
  MoveXYZ(cx - m_center.X(), cy - m_center.Y(), cz - m_center.Z());
}

// Set initial vertex positions (centered at [0, 0, 0])
inline void RectangularPyramid::SetVertex(double wx, double wy, double wz) {
  m_vertices[0].SetXYZ(   0.0,    0.0, 0.0);
  m_vertices[1].SetXYZ(    wx,    0.0, 0.0);
  m_vertices[2].SetXYZ(    wx,     wy, 0.0);
  m_vertices[3].SetXYZ(   0.0,     wy, 0.0);
  m_vertices[4].SetXYZ(0.5*wx, 0.5*wy,  wz);
}

#endif

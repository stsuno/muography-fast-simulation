#ifndef RECTANGULAR_PRISM_H
#define RECTANGULAR_PRISM_H

#include "Simulation/Polyhedron.h"


class RectangularPrism : public Polyhedron {
  public:
    RectangularPrism(double cx, double cy, double cz,
        double wx, double wy, double wz);
    virtual ~RectangularPrism() = default;

  private:
    void SetVertex(double wx, double wy, double wz);
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline RectangularPrism::RectangularPrism(double cx, double cy, double cz,
					  double wx, double wy, double wz) {

  SetFace({
      {0, 1, 2, 3},  // bottom      7------6
      {4, 5, 6, 7},  // top        /|     /|
      {0, 1, 5, 4},  // front     4------5 |   z  y
      {1, 2, 6, 5},  // right     |(3)---|-2   | /
      {2, 3, 7, 6},  // back      |/     |/    |/
      {3, 0, 4, 7}   // left      0------1     o---x
    });

  m_vertices.resize(8);
  SetVertex(wx, wy, wz);
  ComputeCenter();
  MoveXYZ(cx - m_center.X(), cy - m_center.Y(), cz - m_center.Z());
}

// Set initial vertex positions (centered at [0, 0, 0])
inline void RectangularPrism::SetVertex(double wx, double wy, double wz) {
  m_vertices[0].SetXYZ(0.0, 0.0, 0.0);
  m_vertices[1].SetXYZ( wx, 0.0, 0.0);
  m_vertices[2].SetXYZ( wx,  wy, 0.0);
  m_vertices[3].SetXYZ(0.0,  wy, 0.0);
  m_vertices[4].SetXYZ(0.0, 0.0,  wz);
  m_vertices[5].SetXYZ( wx, 0.0,  wz);
  m_vertices[6].SetXYZ( wx,  wy,  wz);
  m_vertices[7].SetXYZ(0.0,  wy,  wz);
}

#endif
